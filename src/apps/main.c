#include "bplcgcd/quantum_channel.h"
#include "bplcgcd/quantum_decoder.h"

#include <errno.h>
#include <inttypes.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

enum { CLI_OK = 0, CLI_USAGE_ERROR = 2, CLI_FAILURE = 4 };

/* Infrequently changed quantum settings stay here. */
static const QuantumLcgcdConfig quantumLcDefaults = {
    .state_limit=4096, .slva_max_candidates=65536, .slva_max_bytes=64*1024*1024,
    .post_policy=BP_POST_ON_FAILURE, .stop_policy=QUANTUM_STOP_TRIVIAL
};
static const double quantumSingletonMessage = 0;

typedef struct {
    const char *hxFile, *hzFile, *rsltFile;
    double pStart, pStop, pStep, nmsScale;
    uint64_t maxFrmNum, stopErrFrmNum, maxIterNum, delta, candidateLimit, randSeed;
} QuantumSimPara;

typedef struct {
    uint64_t errFrmNum, iterNum, searchNum, bpConvFrmNum;
} QuantumStatistic;

static void print_help(FILE *stream)
{
    fputs("BP_LCGCD - Quantum QNMS / BP-LCGCD simulation\n"
          "Usage: bp_lcgcd --quantum HX HZ P_START P_STOP P_STEP MAX_FRAMES MAX_FAILURES\n"
          "       BP_ITER NMS_SCALE DELTA CANDIDATE_LIMIT SEED OUTPUT.csv\n"
          "Or edit run_quantum_sim.sh and run ./run_quantum_sim.sh.\n"
          "\n"
          "Require 0<p_start<=p_stop<1, step>0, frames/iterations>0, 0<scale<=1.\n"
          "MAX_FAILURES=0 disables the failure target; CANDIDATE_LIMIT=0 is unlimited.\n"
          "Both decoders share each physical error and measured syndrome.\n"
          "Exit codes: 0=completed/help, 2=usage error, 4=input/decode/I/O error.\n",stream);
}

static int parse_real(const char *text, double *value)
{
    char *end;
    errno=0;
    *value=strtod(text,&end);
    return end!=text && !*end && errno!=ERANGE && isfinite(*value);
}

static const char *status_name(DecodeStatus status)
{
    static const char *const names[] = {"NOT_RUN", "SUCCESS", "TRUNCATED", "NOT_IMPLEMENTED",
        "INVALID_ARGUMENT", "DIMENSION_MISMATCH", "INVALID_CANDIDATE", "NUMERIC_ERROR",
        "ALLOCATION_FAILED", "INCONSISTENT_SYNDROME"};
    return names[status];
}

static int parse_count(const char *text, uint64_t *value)
{
    char *end;
    errno = 0;
    uintmax_t parsed = strtoumax(text, &end, 10);
    if (*text < '0' || *text > '9' || *end || errno == ERANGE || parsed > UINT64_MAX) return 0;
    *value = (uint64_t)parsed;
    return 1;
}

static int QuantumParaGet(char **argv, QuantumSimPara *para)
{
    *para=(QuantumSimPara){.hxFile=argv[2], .hzFile=argv[3], .rsltFile=argv[14]};
    return parse_real(argv[4],&para->pStart) &&
           parse_real(argv[5],&para->pStop) &&
           parse_real(argv[6],&para->pStep) &&
           parse_count(argv[7],&para->maxFrmNum) &&
           parse_count(argv[8],&para->stopErrFrmNum) &&
           parse_count(argv[9],&para->maxIterNum) &&
           parse_real(argv[10],&para->nmsScale) &&
           parse_count(argv[11],&para->delta) &&
           parse_count(argv[12],&para->candidateLimit) &&
           parse_count(argv[13],&para->randSeed) &&
           para->pStart>0 && para->pStop<1 && para->pStop>=para->pStart &&
           para->pStep>0 && para->pStart+para->pStep>para->pStart &&
           para->maxFrmNum>0 && para->maxIterNum>0 &&
           para->nmsScale>0 && para->nmsScale<=1 &&
           strcmp(para->rsltFile,para->hxFile) && strcmp(para->rsltFile,para->hzFile);
}

static void QuantumResultPrint(FILE *fp, double p, const char *decoder,
                                uint64_t frmNum, const QuantumStatistic *stat)
{
    fprintf(fp,"%.17g,%s,%" PRIu64 ",%" PRIu64 ",%.17g,%" PRIu64 ",%.17g,%" PRIu64 ",%.17g,%" PRIu64 "\n",
        p,decoder,frmNum,stat->errFrmNum,(double)stat->errFrmNum/frmNum,
        stat->iterNum,(double)stat->iterNum/frmNum,
        stat->searchNum,(double)stat->searchNum/frmNum,stat->bpConvFrmNum);
}

int main(int argc, char *argv[])
{
    if (argc==1 || (argc==2 && (!strcmp(argv[1],"--help") || !strcmp(argv[1],"-h")))) {
        print_help(stdout);
        return CLI_OK;
    }
    QuantumSimPara simPara;
    if (argc!=15 || strcmp(argv[1],"--quantum") || !QuantumParaGet(argv,&simPara)) {
        print_help(stderr);
        return CLI_USAGE_ERROR;
    }

    BinaryMatrix hx={0},hz={0};
    QuantumCode quantumCode={0};
    QbpWorkspace qnmsDataMem={0};
    QuantumDecoderWorkspace lcDataMem={0};
    BitVector syndrome={0};
    Pauli *physicalError=NULL;
    double *priorLlr=NULL;
    FILE *hxFp=fopen(simPara.hxFile,"r"), *hzFp=fopen(simPara.hzFile,"r"), *rsltFp=NULL;
    int exitCode=CLI_FAILURE;

    if (!hxFp || !hzFp || binary_matrix_load(hxFp,&hx)!=GF2_OK ||
        binary_matrix_load(hzFp,&hz)!=GF2_OK || quantum_code_init(&quantumCode,&hx,&hz)!=QUANTUM_OK) {
        fputs("INPUT_ERROR: cannot load commuting Hx/Hz matrices.\n",stderr);
        goto finish;
    }
    fclose(hxFp); hxFp=NULL;
    fclose(hzFp); hzFp=NULL;
    size_t n=quantumCode.graph.cols;
    physicalError=malloc(n*sizeof(*physicalError));
    priorLlr=malloc(3*n*sizeof(*priorLlr));
    if (!physicalError || !priorLlr || bit_vector_init(&syndrome,quantumCode.graph.rows)!=GF2_OK ||
        qbp_workspace_init(&qnmsDataMem,&quantumCode)!=QUANTUM_OK) {
        fputs("INPUT_ERROR: quantum workspace allocation failed.\n",stderr);
        goto finish;
    }
    QbpConfig bpPara={.max_iterations=(size_t)simPara.maxIterNum,
        .nms_scale=simPara.nmsScale, .message_clip=quantumSingletonMessage};
    QuantumLcgcdConfig lcPara=quantumLcDefaults;
    lcPara.delta=(size_t)simPara.delta;
    lcPara.candidate_limit=simPara.candidateLimit;
    DecodeStatus status=quantum_decoder_init(&lcDataMem,&quantumCode,&lcPara);
    if (status!=DECODE_SUCCESS) {
        fprintf(stderr,"INPUT_ERROR: LC workspace init %s (check delta and top-of-file defaults).\n",status_name(status));
        goto finish;
    }
    rsltFp=fopen(simPara.rsltFile,"w");
    if (!rsltFp) { perror(simPara.rsltFile); goto finish; }
    fputs("p,decoder,frames,logical_failures,ler,bp_iterations_total,bp_iterations_avg,"
          "lc_searches_total,lc_searches_avg,bp_converged_frames\n",rsltFp);
    QbpPrior prior={n,priorLlr,priorLlr+n,priorLlr+2*n};
    fprintf(stderr,"Quantum: n=%zu m=%zu k=%zu BP_iter=%" PRIu64 " scale=%.6g delta=%" PRIu64
        " candidate_limit=%" PRIu64 " seed=%" PRIu64 "\n",n,quantumCode.graph.rows,
        n-lcDataMem.syndrome_form.rank,simPara.maxIterNum,simPara.nmsScale,
        simPara.delta,simPara.candidateLimit,simPara.randSeed);

    for (size_t pointIdx=0; ; ++pointIdx) {
        double p=simPara.pStart+pointIdx*simPara.pStep;
        if (p>simPara.pStop+simPara.pStep*1e-9) break;
        p=fmin(p,simPara.pStop);
        PauliProbabilities chnnlPara;
        double priorXYZ[3];
        if (quantum_depolarizing_probabilities(p,&chnnlPara)!=QUANTUM_OK ||
            quantum_prior_llr(&chnnlPara,priorXYZ)!=QUANTUM_OK) {
            fputs("INPUT_ERROR: cannot form finite channel prior.\n",stderr);
            goto finish;
        }
        for (size_t i=0; i<n; ++i) {
            priorLlr[i]=priorXYZ[0]; priorLlr[n+i]=priorXYZ[1]; priorLlr[2*n+i]=priorXYZ[2];
        }
        QuantumStatistic qnmsStat={0},lcStat={0};
        uint64_t frmNum=0;
        while (frmNum<simPara.maxFrmNum && (!simPara.stopErrFrmNum ||
               qnmsStat.errFrmNum<simPara.stopErrFrmNum || lcStat.errFrmNum<simPara.stopErrFrmNum)) {
            quantum_depolarizing_sample(p,simPara.randSeed,frmNum,n,physicalError);
            quantum_code_syndrome(&quantumCode,physicalError,n,&syndrome);

            QbpResult qnmsResult;
            QuantumDecodeResult lcResult;
            if (qbp_decode(&qnmsDataMem,&bpPara,&prior,&syndrome,&qnmsResult)!=QUANTUM_OK) {
                fprintf(stderr,"DECODE_ERROR: QNMS p=%.6g frame=%" PRIu64 "\n",p,frmNum);
                goto finish;
            }
            status=quantum_bp_lcgcd_decode(&lcDataMem,&bpPara,&prior,&syndrome,&lcResult);
            if (status!=DECODE_SUCCESS && status!=DECODE_TRUNCATED) {
                fprintf(stderr,"DECODE_ERROR: BP-LCGCD %s\n",status_name(status));
                goto finish;
            }
            QuantumRecoveryCheck qnmsCheck,lcCheck={0};
            quantum_recovery_check(&lcDataMem,physicalError,qnmsResult.hard_estimate,n,&qnmsCheck);
            if (lcResult.candidate_valid)
                quantum_recovery_check(&lcDataMem,physicalError,lcResult.recovery,n,&lcCheck);

            ++frmNum;
            qnmsStat.errFrmNum+=!qnmsCheck.logical_success;
            qnmsStat.iterNum+=qnmsResult.iterations_used;
            qnmsStat.bpConvFrmNum+=qnmsResult.syndrome_matched;
            lcStat.errFrmNum+=!lcCheck.logical_success;
            lcStat.iterNum+=lcResult.posterior.iterations_used;
            lcStat.searchNum+=lcResult.stats.candidates_evaluated;
            lcStat.bpConvFrmNum+=lcResult.posterior.syndrome_matched;
        }
        fprintf(stderr,"p=%.6g frames=%" PRIu64 " | QNMS fail=%" PRIu64 " LER=%.6g BP=%.3f"
            " | BP-LCGCD fail=%" PRIu64 " LER=%.6g BP=%.3f LC=%.3f\n",p,frmNum,
            qnmsStat.errFrmNum,(double)qnmsStat.errFrmNum/frmNum,(double)qnmsStat.iterNum/frmNum,
            lcStat.errFrmNum,(double)lcStat.errFrmNum/frmNum,(double)lcStat.iterNum/frmNum,(double)lcStat.searchNum/frmNum);
        QuantumResultPrint(rsltFp,p,"QNMS",frmNum,&qnmsStat);
        QuantumResultPrint(rsltFp,p,"BP-LCGCD",frmNum,&lcStat);
        if (ferror(rsltFp) || fflush(rsltFp)!=0) { fputs("IO_ERROR: CSV write failed.\n",stderr); goto finish; }
    }
    exitCode=CLI_OK;
finish:
    if (hxFp) fclose(hxFp);
    if (hzFp) fclose(hzFp);
    if (rsltFp && fclose(rsltFp)!=0) exitCode=CLI_FAILURE;
    quantum_decoder_free(&lcDataMem);
    qbp_workspace_free(&qnmsDataMem);
    quantum_code_free(&quantumCode);
    binary_matrix_free(&hx); binary_matrix_free(&hz);
    bit_vector_free(&syndrome);
    free(physicalError); free(priorLlr);
    return exitCode;
}
