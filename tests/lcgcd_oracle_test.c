#include "bplcgcd/decoder.h"

#include <float.h>
#include <inttypes.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

enum { MAX_N=10, MAX_ROWS=6 };
static unsigned checks, fixtures, cases, calls, proven, truncated, inconsistent, replays;
static unsigned exact_proven, exact_truncated, exact_inconsistent;
static unsigned delta_zero, delta_full, delta_middle, decimal_cases, rounded_scores, numeric_ml_mismatches;
static long double max_score_error, max_ml_excess;
static long double max_relative_excess;
static uint32_t seed=20261124u;
static const char *debug_path;
static const DecodeInput *current_input;
static size_t current_delta;
static uint64_t current_limit;

static void failure(const char *condition, int line)
{
    fprintf(stderr,"S15 failure line=%d fixture=%u delta=%zu limit=%" PRIu64 ": %s\n",
            line,fixtures,current_delta,current_limit,condition);
    if (current_input) {
        const DecodeInput *in=current_input;
        fprintf(stderr,"mode=%d H=%zux%zu, seed=20261124\n",(int)in->mode,in->h->rows,in->length);
        for (size_t i=0; i<in->h->rows; ++i) {
            for (size_t j=0; j<in->length; ++j) fprintf(stderr,"%u ",in->h->data[i*in->length+j]);
            fputc('\n',stderr);
        }
        fputs("soft:",stderr);
        for (size_t j=0; j<in->length; ++j) fprintf(stderr," %.17g",in->llr ? in->llr[j] : in->weights[j]);
        fputs("\ntarget:",stderr);
        if (in->s_target) for (size_t i=0; i<in->s_target->length; ++i) fprintf(stderr,"%u",in->s_target->data[i]);
        fputs(" hint:",stderr);
        if (in->z_hint) for (size_t j=0; j<in->length; ++j) fprintf(stderr,"%u",in->z_hint->data[j]);
        fputc('\n',stderr);
    }
    exit(EXIT_FAILURE);
}
#define CHECK(condition) do { if (!(condition)) failure(#condition,__LINE__); ++checks; } while (0)

static unsigned random_value(void)
{
    seed=1664525u*seed+1013904223u;
    return seed>>16;
}

/* Original H and vector parity only; no GF(2), S03, S10 or search helpers. */
static unsigned syndrome(const BinaryMatrix *h, unsigned word)
{
    unsigned result=0;
    for (size_t i=0; i<h->rows; ++i) {
        unsigned parity=0;
        for (size_t j=0; j<h->cols; ++j) parity ^= h->data[i*h->cols+j] & ((word>>j)&1u);
        result |= parity<<i;
    }
    return result;
}

static unsigned pack(const BitVector *v)
{
    unsigned word=0;
    for (size_t j=0; j<v->length; ++j) { CHECK(v->data[j]<=1); word |= (unsigned)v->data[j]<<j; }
    return word;
}

static long double weight(unsigned f, const double *weights, size_t n)
{
    long double cost=0;
    for (size_t j=0; j<n; ++j) if ((f>>j)&1u) cost+=(long double)weights[j];
    return cost;
}

static long double tolerance(long double a, long double b, size_t n, int exact)
{
    /* Strict families use zero tolerance. Decimal diagnostics use a stated
     * forward-error bound and separately report higher-precision ML misses. */
    long double scale=fabsl(a)>fabsl(b) ? fabsl(a) : fabsl(b);
    return exact ? 0 : 4*(n+1)*DBL_EPSILON*scale;
}

static DecodeStatus run(const DecodeConfig *config, const DecodeInput *input, DecodeResult *result)
{
    current_input=input; current_delta=config->delta; current_limit=config->candidate_limit;
    /* Optional per-invocation log keeps user debug prints and bounds test I/O.
     * Only this test's stdout is redirected; summaries/failures use stderr. */
    if (debug_path) CHECK(freopen(debug_path,"w",stdout));
    printf("S15 fixture=%u engine=%s delta=%zu limit=%" PRIu64 "\n",
           fixtures,config->decoder_name,config->delta,config->candidate_limit);
    ++calls;
    return decoder_decode(config,input,result);
}

static void same_result(const DecodeResult *a, const DecodeResult *b)
{
    CHECK(a->status==b->status && a->candidate_valid==b->candidate_valid);
    CHECK(a->optimality_proven==b->optimality_proven && a->stop_reason==b->stop_reason);
    if (a->candidate_valid) CHECK(a->cost==b->cost && !memcmp(a->estimate.data,b->estimate.data,a->estimate.length));
    CHECK(a->stats.candidates_popped==b->stats.candidates_popped);
    CHECK(a->stats.candidates_evaluated==b->stats.candidates_evaluated);
    CHECK(a->stats.nodes_expanded==b->stats.nodes_expanded && a->stats.batches==b->stats.batches);
}

static void verify_input(DecodeInput *input, int exact)
{
    const BinaryMatrix *h=input->h;
    unsigned reference=0, target=input->mode==DECODE_CLASSICAL ? 0 : pack(input->s_target);
    double weights[MAX_N];
    for (size_t j=0; j<h->cols; ++j) {
        unsigned bit=input->mode==DECODE_SYNDROME ? (input->z_hint ? input->z_hint->data[j] : 0) : input->llr[j]<0;
        reference |= bit<<j;
        weights[j]=input->weights ? input->weights[j] : fabs(input->llr[j]);
    }
    long double optimum=INFINITY;
    unsigned minimizers=0;
    for (unsigned f=0; f<(1u<<h->cols); ++f) if (syndrome(h,reference^f)==target) {
        long double cost=weight(f,weights,h->cols);
        if (cost<optimum) { optimum=cost; minimizers=1; }
        else if (cost==optimum) ++minimizers;
    }
    unsigned total=1u<<(h->cols-input->systematic->rank);
    uint64_t limits[4]={1,2,total,0};
    for (size_t delta=0; delta<=input->systematic->rank; ++delta) {
        long double previous=INFINITY;
        for (size_t trial=0; trial<4; ++trial) {
            uint64_t limit=limits[trial];
            if ((trial==1 && total<2) || (trial==2 && total<=2)) continue;
            uint8_t bits[MAX_N], other_bits[MAX_N];
            DecodeResult result={.estimate={h->cols,bits}}, other={.estimate={h->cols,other_bits}};
            DecodeConfig config={.decoder_name="lcgcd",.delta=delta,.candidate_limit=limit};
            DecodeStatus status=run(&config,input,&result);
            ++cases;
            if (!delta) ++delta_zero;
            else if (delta==input->systematic->rank) ++delta_full;
            else ++delta_middle;
            CHECK(status==result.status && result.stats.batches==0);
            CHECK(result.stats.candidates_popped==result.stats.candidates_evaluated);
            CHECK(!limit || result.stats.candidates_evaluated<=limit);
            CHECK(result.stats.nodes_expanded<=result.stats.candidates_popped);
            if (!isfinite(optimum)) {
                ++inconsistent;
                if (exact) ++exact_inconsistent;
                CHECK(status==DECODE_INCONSISTENT_SYNDROME && !result.candidate_valid && !result.optimality_proven);
                CHECK(result.stop_reason==DECODE_STOP_NONE && !result.stats.candidates_popped);
            } else {
                CHECK(status==DECODE_SUCCESS || status==DECODE_TRUNCATED);
                CHECK(result.candidate_valid && syndrome(h,pack(&result.estimate))==target);
                long double actual=weight(pack(&result.estimate)^reference,weights,h->cols);
                long double error=fabsl(actual-result.cost);
                CHECK(error<=tolerance(actual,result.cost,h->cols,exact));
                CHECK(actual>=optimum);
                CHECK(!isfinite(previous) || actual<=previous+tolerance(actual,previous,h->cols,exact));
                previous=actual;
                if (!exact) {
                    ++decimal_cases; rounded_scores += error!=0;
                    if (error>max_score_error) max_score_error=error;
                }
                if (status==DECODE_SUCCESS) {
                    ++proven;
                    if (exact) ++exact_proven;
                    CHECK(result.optimality_proven);
                    CHECK(result.stop_reason==DECODE_STOP_BOUND || result.stop_reason==DECODE_STOP_EXHAUSTED);
                    CHECK(actual-optimum<=tolerance(actual,optimum,h->cols,exact));
                    if (!exact && actual-optimum>max_ml_excess) max_ml_excess=actual-optimum;
                    if (!exact && actual>optimum) ++numeric_ml_mismatches;
                    if (!exact && actual>0 && (actual-optimum)/actual>max_relative_excess)
                        max_relative_excess=(actual-optimum)/actual;
                } else {
                    ++truncated;
                    if (exact) ++exact_truncated;
                    CHECK(limit && result.stats.candidates_evaluated==limit);
                    CHECK(!result.optimality_proven && result.stop_reason==DECODE_STOP_BUDGET);
                }
                if (!limit || limit>=total || delta==input->systematic->rank) CHECK(status==DECODE_SUCCESS);
                if (delta==input->systematic->rank) {
                    CHECK(result.stats.candidates_evaluated==1);
                    if (delta) CHECK(!result.stats.nodes_expanded);
                }
            }
            if (!delta) {
                config.decoder_name="gcd";
                run(&config,input,&other);
                same_result(&result,&other); /* Exact delta=0 reduction, even at ties/limits. */
            }
            if (minimizers>1 && (limit==0 || limit==1)) {
                config.decoder_name="lcgcd";
                run(&config,input,&other);
                same_result(&result,&other);
                ++replays;
            }
        }
    }
}

static void exercise(uint8_t *data, size_t rows, size_t n, const double *weights,
                     const double *llr, unsigned hint_word, unsigned target_word, int exact)
{
    ++fixtures;
    BinaryMatrix h={rows,n,data};
    BinarySystematic form={0};
    uint8_t original[MAX_ROWS*MAX_N], reduced[MAX_ROWS*MAX_N], transform[MAX_ROWS*MAX_ROWS];
    size_t columns[MAX_N], inverse[MAX_N];
    double saved_weights[MAX_N], saved_llr[MAX_N];
    memcpy(original,data,rows*n);
    memcpy(saved_weights,weights,n*sizeof(*weights));
    memcpy(saved_llr,llr,n*sizeof(*llr));
    if (debug_path) CHECK(freopen(debug_path,"w",stdout));
    CHECK(binary_matrix_systematize(&h,&form)==GF2_OK);
    /* Independent rank: size of the image of H, enumerating all original words. */
    unsigned images[1u<<MAX_ROWS]={0}, image_count=0;
    for (unsigned f=0; f<(1u<<n); ++f) images[syndrome(&h,f)]=1;
    for (unsigned s=0; s<(1u<<rows); ++s) image_count+=images[s];
    CHECK(image_count==(1u<<form.rank));
    memcpy(reduced,form.reduced.data,rows*n);
    memcpy(transform,form.row_transform.data,rows*rows);
    memcpy(columns,form.sys_to_orig,n*sizeof(*columns));
    memcpy(inverse,form.orig_to_sys,n*sizeof(*inverse));
    uint8_t target_bits[MAX_ROWS], hint_bits[MAX_N];
    for (size_t i=0; i<rows; ++i) target_bits[i]=(uint8_t)((target_word>>i)&1u);
    for (size_t j=0; j<n; ++j) hint_bits[j]=(uint8_t)((hint_word>>j)&1u);
    BitVector target={rows,target_bits}, hint={n,hint_bits};
    for (unsigned mode=0; mode<4; ++mode) {
        DecodeInput input={.mode=mode==0 ? DECODE_CLASSICAL : mode==3 ? DECODE_BP_SYNDROME : DECODE_SYNDROME,
            .h=&h,.length=n,.systematic=&form};
        if (input.mode==DECODE_SYNDROME) { input.weights=weights; if (mode==2) input.z_hint=&hint; }
        else input.llr=llr;
        if (mode) input.s_target=&target;
        current_input=&input;
        verify_input(&input,exact);
    }
    current_input=NULL;
    CHECK(!memcmp(original,data,rows*n));
    CHECK(!memcmp(reduced,form.reduced.data,rows*n));
    CHECK(!memcmp(transform,form.row_transform.data,rows*rows));
    CHECK(!memcmp(columns,form.sys_to_orig,n*sizeof(*columns)));
    CHECK(!memcmp(inverse,form.orig_to_sys,n*sizeof(*inverse)));
    CHECK(!memcmp(saved_weights,weights,n*sizeof(*weights)) && !memcmp(saved_llr,llr,n*sizeof(*llr)));
    CHECK(pack(&target)==(target_word&((1u<<rows)-1)) && pack(&hint)==(hint_word&((1u<<n)-1)));
    binary_systematic_free(&form);
}

static void delta_zero_regression(void)
{
    uint8_t data[]={1,0,1,1, 0,1,0,1};
    double weights[]={10,10,1,1}, llr[]={10,-10,1,-1};
    exercise(data,2,4,weights,llr,5,1,1);
}

int main(int argc, char **argv)
{
    if (argc>1) debug_path=argv[1];
    delta_zero_regression();
    /* All 64 binary 2x3 H and every external syndrome; includes impossible targets. */
    for (unsigned matrix=0; matrix<64; ++matrix) for (unsigned target=0; target<4; ++target) {
        uint8_t data[6]; double weights[3], llr[3];
        for (unsigned j=0; j<6; ++j) data[j]=(uint8_t)((matrix>>j)&1u);
        for (unsigned j=0; j<3; ++j) {
            weights[j]=matrix%4==0 ? 0 : matrix%4==1 ? 1 : matrix%4==2 ? (double)((3*j+matrix)%7)/4 : (double)(1u<<j);
            llr[j]=(matrix>>j)&1u ? -weights[j] : weights[j];
        }
        exercise(data,2,3,weights,llr,5,target,1);
    }
    for (unsigned trial=0; trial<30; ++trial) {
        size_t n=1+trial%MAX_N, rows=1+trial%4;
        uint8_t data[MAX_ROWS*MAX_N]={0}, permuted[MAX_ROWS*MAX_N], redundant[MAX_ROWS*MAX_N];
        double weights[MAX_N], llr[MAX_N], pw[MAX_N], pl[MAX_N];
        size_t permutation[MAX_N];
        for (size_t j=0; j<rows*n; ++j) data[j]=(uint8_t)(random_value()&1u);
        for (size_t j=0; j<n; ++j) {
            weights[j]=trial%5==0 ? 0 : trial%5==1 ? 1 : (double)(random_value()%25)/4;
            llr[j]=random_value()&1u ? -weights[j] : weights[j];
            permutation[j]=j;
        }
        for (size_t j=n; j>1; --j) {
            size_t p=random_value()%j, saved=permutation[j-1]; permutation[j-1]=permutation[p]; permutation[p]=saved;
        }
        unsigned word=random_value()&((1u<<n)-1), hint=(random_value()|1u)&((1u<<n)-1);
        BinaryMatrix h={rows,n,data};
        unsigned target=syndrome(&h,word), permuted_hint=0;
        exercise(data,rows,n,weights,llr,hint,target,1);
        for (size_t j=0; j<n; ++j) {
            pw[j]=weights[permutation[j]]; pl[j]=llr[permutation[j]];
            permuted_hint |= ((hint>>permutation[j])&1u)<<j;
            for (size_t i=0; i<rows; ++i) permuted[i*n+j]=data[i*n+permutation[j]];
        }
        exercise(permuted,rows,n,pw,pl,permuted_hint,target,1);
        memcpy(redundant,data,rows*n);
        for (size_t j=0; j<n; ++j) { redundant[rows*n+j]=data[j]^data[(rows-1)*n+j]; redundant[(rows+1)*n+j]=0; }
        unsigned extra=(target&1u)^((target>>(rows-1))&1u);
        unsigned redundant_target=target|(extra<<rows);
        if (trial%3==0) redundant_target |= 1u<<(rows+1); /* Impossible zero-row RHS. */
        exercise(redundant,rows+2,n,weights,llr,hint,redundant_target,1);
    }
    uint8_t identity[]={1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
    const double binary[]={1,2,4,8}, signed_binary[]={-1,2,-4,8};
    exercise(identity,4,4,binary,signed_binary,5,15,1); /* K=0, including L=0. */
    uint8_t hamming[28];
    for (size_t i=0; i<3; ++i) for (size_t j=0; j<7; ++j) hamming[i*7+j]=(uint8_t)(((j+1)>>i)&1u);
    for (size_t j=0; j<7; ++j) hamming[21+j]=hamming[j]^hamming[7+j];
    const double hw[]={8,3,6,1,2,4,7}, hl[]={-8,-3,-6,1,2,4,7};
    exercise(hamming,4,7,hw,hl,37,7,1); /* Self-authored Hamming [7,4], redundant row. */
    const double decimal[]={0.1,0.2,0.3,0.4}, signed_decimal[]={-0.1,0.2,-0.3,0.4};
    exercise(identity,4,4,decimal,signed_decimal,5,15,0);
    const double dynamic[]={1e15,0.1,0.2,0.3}, signed_dynamic[]={-1e15,0.1,-0.2,0.3};
    exercise(identity,4,4,dynamic,signed_dynamic,5,15,0);
    uint8_t free_case[]={1,0,1,0,1, 0,1,0,1,1};
    const double fw[]={0.1,0.2,0.3,0.4,0.3}, fl[]={-0.1,0.2,-0.3,0.4,0.3};
    exercise(free_case,2,5,fw,fl,21,3,0);
    uint8_t near_tie[]={1,0,0, 0,1,1};
    const double nw[]={1e16,0.2,0.1}, nl[]={-1e16,0.2,0.1};
    exercise(near_tie,2,3,nw,nl,0,3,0); /* Tiny score differences lost after a mandatory large term. */
    CHECK(delta_zero && delta_full && delta_middle && proven && truncated && inconsistent && replays && rounded_scores);
    fprintf(stderr,"S15 validation: checks=%u fixtures=%u LC_cases=%u decoder_calls=%u\n",checks,fixtures,cases,calls);
    fprintf(stderr,"exact_cases=%u exact_proven=%u exact_truncated=%u exact_inconsistent=%u\n",
            cases-decimal_cases,exact_proven,exact_truncated,exact_inconsistent);
    fprintf(stderr,"delta_zero=%u delta_full=%u delta_middle=%u proven=%u truncated=%u inconsistent=%u tie_replays=%u\n",
            delta_zero,delta_full,delta_middle,proven,truncated,inconsistent,replays);
    fprintf(stderr,"decimal_cases=%u rounded_scores=%u numeric_ML_mismatches=%u max_score_error=%.9Lg max_ML_excess=%.9Lg relative_excess=%.9Lg LDBL_MANT_DIG=%d\n",
            decimal_cases,rounded_scores,numeric_ml_mismatches,max_score_error,max_ml_excess,max_relative_excess,LDBL_MANT_DIG);
    return EXIT_SUCCESS;
}
