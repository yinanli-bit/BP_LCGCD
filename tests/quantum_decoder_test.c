#include "bplcgcd/quantum_decoder.h"
#include "bplcgcd/quantum_channel.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

static unsigned checks, oracle_cases, logical_cases;
static unsigned bound_stops, bound_saved;
static const char *debug_log;
#define CHECK(c) do { if (!(c)) { \
    fprintf(stderr,"%s:%d: %s\n",__FILE__,__LINE__,#c); exit(EXIT_FAILURE); \
} ++checks; } while (0)

#ifdef QBP_TRACK_ALLOCATIONS
static int tracking;
static unsigned allocations, frees, eliminations;
void *__real_malloc(size_t);
void *__real_calloc(size_t,size_t);
void *__real_realloc(void *,size_t);
void __real_free(void *);
Gf2Status __real_binary_matrix_systematize(const BinaryMatrix *,BinarySystematic *);
void *__wrap_malloc(size_t n) { if (tracking) ++allocations; return __real_malloc(n); }
void *__wrap_calloc(size_t n,size_t s) { if (tracking) ++allocations; return __real_calloc(n,s); }
void *__wrap_realloc(void *p,size_t n) { if (tracking) ++allocations; return __real_realloc(p,n); }
void __wrap_free(void *p) { if (tracking) ++frees; __real_free(p); }
Gf2Status __wrap_binary_matrix_systematize(const BinaryMatrix *m,BinarySystematic *s)
{
    if (tracking) ++eliminations;
    return __real_binary_matrix_systematize(m,s);
}
static void track(int enabled) { tracking=enabled; }
#else
static void track(int enabled) { (void)enabled; }
#endif

enum { MAX_N=5, MAX_M=4, PATTERNS=1024 };
static const char labels[] = "IXYZ";
static const unsigned anti[4][4] = {{0,0,0,0},{0,0,1,1},{0,1,0,1},{0,1,1,0}};
typedef struct { size_t m,n; const char *rows; } Fixture;
static const Fixture toys[] = {
    {1,1,"Z"}, {3,2,"XZ" "ZX" "YY"}, {2,3,"XYZ" "ZXI"},
    {3,3,"ZZI" "IZZ" "ZIZ"}, {1,3,"III"}
};
static const Fixture five = {4,5,"XZZXI" "ZZXIX" "ZXIXZ" "XIXZZ"};

static int near(double a,double b) { return fabs(a-b) <= 2e-10*fmax(1,fabs(b)); }
static Pauli label(char c) { return (Pauli)(strchr(labels,c)-labels); }
static unsigned pack(const Pauli *p,size_t n)
{
    unsigned word=0;
    for (size_t i=0;i<n;++i) word |= (unsigned)p[i] << (2*i);
    return word;
}
static void unpack(unsigned word,Pauli *p,size_t n)
{
    for (size_t i=0;i<n;++i) p[i]=(Pauli)((word>>(2*i))&3u);
}
static void reset_log(void)
{
    if (debug_log) CHECK(freopen(debug_log,"w",stdout) != NULL);
}

/* Integer Pauli/group oracle: no GF2/systematic/local/SLVA production helper.
 * With labels I=0,X=1,Y=2,Z=3, XOR is Pauli multiplication modulo phase. */
static unsigned oracle_syndrome(const Fixture *f,unsigned word)
{
    unsigned syndrome=0;
    for (size_t r=0;r<f->m;++r) {
        unsigned parity=0;
        for (size_t i=0;i<f->n;++i) parity ^= anti[(word>>(2*i))&3u][label(f->rows[r*f->n+i])];
        syndrome |= parity << r;
    }
    return syndrome;
}
static size_t oracle_group(const Fixture *f,uint8_t *member)
{
    memset(member,0,PATTERNS);
    for (unsigned subset=0;subset<(1u<<f->m);++subset) {
        unsigned word=0;
        for (size_t r=0;r<f->m;++r)
            if ((subset>>r)&1u)
                for (size_t i=0;i<f->n;++i) word ^= (unsigned)label(f->rows[r*f->n+i]) << (2*i);
        member[word]=1;
    }
    size_t count=0;
    for (unsigned word=0;word<(1u<<(2*f->n));++word) count += member[word];
    return count;
}
static double oracle_score(unsigned word,const QbpResult *post)
{
    long double total=0;
    for (size_t i=0;i<post->length;++i) {
        long double mass[] = {1,expl(-post->post_x[i]),expl(-post->post_y[i]),expl(-post->post_z[i])};
        long double most=mass[0];
        for (size_t p=1;p<4;++p) if (mass[p]>most) most=mass[p];
        total += logl(most/mass[(word>>(2*i))&3u]);
    }
    return (double)total;
}
static size_t oracle_best(const Fixture *f,const QbpResult *post,unsigned syndrome,double *best)
{
    size_t count=0;
    *best=INFINITY;
    for (unsigned word=0;word<(1u<<(2*f->n));++word) {
        if (oracle_syndrome(f,word) != syndrome) continue;
        double cost=oracle_score(word,post);
        if (cost<*best) *best=cost;
        ++count;
    }
    return count;
}

static QuantumCode make_code(const Fixture *f)
{
    uint8_t x[MAX_M*MAX_N],z[MAX_M*MAX_N];
    CHECK(f->m<=MAX_M && f->n<=MAX_N);
    for (size_t k=0;k<f->m*f->n;++k) {
        x[k]=f->rows[k]=='X' || f->rows[k]=='Y';
        z[k]=f->rows[k]=='Z' || f->rows[k]=='Y';
    }
    BinaryMatrix hx={f->m,f->n,x},hz={f->m,f->n,z};
    QuantumCode code={0};
    CHECK(quantum_code_init(&code,&hx,&hz)==QUANTUM_OK);
    return code;
}
static QuantumLcgcdConfig config(size_t delta)
{
    return (QuantumLcgcdConfig){.delta=delta,.state_limit=32,
        .slva_max_candidates=512,.slva_max_bytes=1024*1024,.stop_policy=QUANTUM_STOP_DISABLED};
}
static QbpResult posterior(double values[3][MAX_N],size_t n)
{
    return (QbpResult){.length=n,.post_x=values[0],.post_y=values[1],.post_z=values[2]};
}
static DecodeStatus run_post(QuantumDecoderWorkspace *w,const QbpResult *post,
                             const BitVector *syn,QuantumDecodeResult *out)
{
    track(1);
    DecodeStatus s=quantum_lcgcd_decode(w,post,syn,out);
    track(0);
    return s;
}
static DecodeStatus run_bp(QuantumDecoderWorkspace *w,const QbpConfig *bp,const QbpPrior *prior,
                           const BitVector *syn,QuantumDecodeResult *out)
{
    track(1);
    DecodeStatus s=quantum_bp_lcgcd_decode(w,bp,prior,syn,out);
    track(0);
    return s;
}

/* Read the remaining SLVA groups only to identify their domains. Syndrome and
 * Pauli scores are still computed independently, without a search helper. */
static int in_group(unsigned word,const QuantumDecoderWorkspace *w,const SlvaCandidate *group)
{
    size_t n=w->code->graph.cols;
    for (size_t j=group->free_length;j<w->trellis.length;++j) {
        size_t col=w->trellis.original_columns[j],i=col%n;
        unsigned p=(word>>(2*i))&3u;
        unsigned bit=col<n ? (p==PAULI_X || p==PAULI_Y) : (p==PAULI_Z || p==PAULI_Y);
        if ((bit^w->reference.data[col]) != group->bits[j]) return 0;
    }
    return 1;
}
static void verify_pauli_bound(const Fixture *f,const QbpResult *post,unsigned target,
                               const QuantumDecoderWorkspace *w,const QuantumDecodeResult *out)
{
    double remaining=INFINITY;
    for (unsigned word=0;word<(1u<<(2*f->n));++word) {
        double binary=0,full=0;
        for (size_t j=0;j<w->trellis.length;++j) {
            size_t col=w->trellis.original_columns[j],i=col%f->n;
            unsigned p=(word>>(2*i))&3u;
            unsigned bit=col<f->n ? (p==PAULI_X || p==PAULI_Y) : (p==PAULI_Z || p==PAULI_Y);
            if (bit^w->reference.data[col]) binary += w->reliability[col];
        }
        for (size_t i=0;i<f->n;++i) full += w->swt[i].w[(word>>(2*i))&3u];
        CHECK((long double)binary-w->stop_correction <= (long double)full);
        if (oracle_syndrome(f,word) != target) continue;
        int unseen=w->list.pending.bits && in_group(word,w,&w->list.pending);
        for (size_t g=0;!unseen && g<w->list.size;++g) unseen=in_group(word,w,&w->list.heap[g]);
        if (unseen) remaining=fmin(remaining,oracle_score(word,post));
    }
    if (out->early_stop!=QUANTUM_EARLY_NONE)
        CHECK(out->remaining_lower_bound <= remaining+2e-12);
}

static void test_soft_information(void)
{
    double probabilities[3][4]={{0.5,0.05,0.4,0.05},{0.3,0.35,0.05,0.3},{0.25,0.25,0.25,0.25}};
    double values[3][MAX_N]={{0}};
    for (size_t i=0;i<3;++i)
        for (size_t p=0;p<3;++p) values[p][i]=log(probabilities[i][0]/probabilities[i][p+1]);
    QbpResult post=posterior(values,3);
    PauliSwt swt[3]; PartialSwt partial[3];
    uint8_t hard_data[6]; double reliability[6];
    BitVector hard={6,hard_data};
    CHECK(qbp_posterior_to_swt(&post,swt)==QUANTUM_OK);
    CHECK(quantum_swt_to_partial(swt,3,partial,&hard,reliability)==QUANTUM_OK);
    for (size_t i=0;i<3;++i) {
        double most=0;
        unsigned hard_pauli=0;
        for (unsigned p=0;p<4;++p)
            if (probabilities[i][p]>most) { most=probabilities[i][p]; hard_pauli=p; }
        unsigned minimum=0;
        for (unsigned p=0;p<4;++p) {
            CHECK(near(swt[i].w[p],log(most/probabilities[i][p])));
            if (swt[i].w[p]<swt[i].w[minimum]) minimum=p;
            for (unsigned q=0;q<4;++q)
                CHECK(near(swt[i].w[p]-swt[i].w[q],log(probabilities[i][q]/probabilities[i][p])));
        }
        CHECK(minimum==hard_pauli);
        const double *q=probabilities[i];
        CHECK(near(partial[i].x[0],log(most/(q[0]+q[3]))));
        CHECK(near(partial[i].x[1],log(most/(q[1]+q[2]))));
        CHECK(near(partial[i].z[0],log(most/(q[0]+q[1]))));
        CHECK(near(partial[i].z[1],log(most/(q[3]+q[2]))));
        CHECK(near(reliability[i],fabs(log((q[0]+q[3])/(q[1]+q[2])))));
        CHECK(near(reliability[3+i],fabs(log((q[0]+q[1])/(q[3]+q[2])))));
    }
    CHECK(swt[1].w[PAULI_X]==0 && hard_data[1]==0 && hard_data[4]==0); /* Marginal hard != Pauli hard. */
    CHECK(reliability[0]+reliability[3]>swt[0].w[PAULI_Y]); /* Not a valid full-cost lower bound. */
    CHECK(near(swt[0].w[PAULI_Y],log(1.25)));
    CHECK(!near(swt[0].w[PAULI_Y],swt[0].w[PAULI_X]+swt[0].w[PAULI_Z]));
    values[0][0]=-1000; values[1][0]=1000; values[2][0]=0;
    post=posterior(values,1); hard.length=2;
    CHECK(qbp_posterior_to_swt(&post,swt)==QUANTUM_OK);
    CHECK(quantum_swt_to_partial(swt,1,partial,&hard,reliability)==QUANTUM_OK);
    CHECK(swt[0].w[0]==1000 && swt[0].w[1]==0 && swt[0].w[2]==2000);
    CHECK(isfinite(reliability[0]) && isfinite(reliability[1]) && hard_data[0]==1 && hard_data[1]==0);
}

static void test_exhaustive_postprocessing(void)
{
    for (size_t fixture=0;fixture<sizeof(toys)/sizeof(*toys);++fixture) {
        const Fixture *f=&toys[fixture];
        uint8_t members[PATTERNS];
        size_t group_size=oracle_group(f,members),rank=0;
        while (((size_t)1<<rank)<group_size) ++rank;
        QuantumCode code=make_code(f);
        for (size_t delta=0;delta<=rank;++delta) {
            QuantumDecoderWorkspace w={0};
            QuantumLcgcdConfig cfg=config(delta);
            CHECK(quantum_decoder_init(&w,&code,&cfg)==DECODE_SUCCESS);
            QuantumDecoderWorkspace early={0};
            cfg.stop_policy=QUANTUM_STOP_TRIVIAL;
            CHECK(quantum_decoder_init(&early,&code,&cfg)==DECODE_SUCCESS);
            CHECK(w.list.path_pool && w.syndrome_form.rank==rank && w.stabilizer_form.rank==rank);
            void *pool=w.list.path_pool, *cells=w.list.cells, *messages=w.bp.messages;
            for (unsigned trial=0;trial<3;++trial) {
                double values[3][MAX_N]={{0}};
                for (size_t p=0;p<3;++p)
                    for (size_t i=0;i<f->n;++i)
                        values[p][i]=trial ? ((int)((trial*11+p*3+i*7)%17)-8)*0.25 : 0;
                QbpResult post=posterior(values,f->n);
                for (unsigned target=0;target<(1u<<f->m);++target) {
                    reset_log();
                    uint8_t syn_data[MAX_M];
                    for (size_t r=0;r<f->m;++r) syn_data[r]=(target>>r)&1u;
                    BitVector syn={f->m,syn_data};
                    double best;
                    size_t count=oracle_best(f,&post,target,&best);
                    QuantumDecodeResult out;
                    DecodeStatus status=run_post(&w,&post,&syn,&out);
                    if (!count) {
                        CHECK(status==DECODE_INCONSISTENT_SYNDROME);
                        CHECK(!out.candidate_valid && !out.optimality_proven && out.stats.candidates_evaluated==0);
                    } else {
                        CHECK(status==DECODE_SUCCESS && out.status==status);
                        CHECK(out.candidate_valid && out.syndrome_matched && out.optimality_proven);
                        CHECK(out.stop_reason==DECODE_STOP_EXHAUSTED && out.postprocess_invoked);
                        CHECK(oracle_syndrome(f,pack(out.recovery,f->n))==target);
                        CHECK(near(out.cost,best) && near(out.cost,oracle_score(pack(out.recovery,f->n),&post)));
                        CHECK(out.stats.candidates_popped==count && out.stats.candidates_evaluated==count);
                        CHECK(out.stats.partial_patterns_generated==count && out.stats.trellis_states_visited>0);
                    }
                    CHECK(w.list.path_pool==pool && w.list.cells==cells && w.bp.messages==messages);
                    QuantumDecodeResult stopped;
                    CHECK(run_post(&early,&post,&syn,&stopped)==status);
                    if (count) {
                        CHECK(stopped.candidate_valid && stopped.syndrome_matched && stopped.optimality_proven);
                        CHECK(near(stopped.cost,best) && stopped.cost==out.cost);
                        CHECK(memcmp(stopped.recovery,out.recovery,f->n*sizeof(Pauli))==0);
                        CHECK(stopped.stats.candidates_evaluated<=count);
                        if (stopped.early_stop==QUANTUM_EARLY_TRIVIAL) {
                            CHECK(stopped.stop_reason==DECODE_STOP_BOUND && stopped.remaining_lower_bound>=stopped.cost);
                            ++bound_stops;
                            bound_saved += (unsigned)(count-stopped.stats.candidates_evaluated);
                        } else CHECK(stopped.stop_reason==DECODE_STOP_EXHAUSTED);
                        verify_pauli_bound(f,&post,target,&early,&stopped);
                    } else CHECK(!stopped.candidate_valid && !stopped.optimality_proven);
                    ++oracle_cases;
                }
            }
            quantum_decoder_free(&w);
            quantum_decoder_free(&early);
        }
        quantum_code_free(&code);
    }
}

static void test_stopping_policies(void)
{
    const Fixture f={2,2,"XX" "ZZ"};
    QuantumCode code=make_code(&f);
    uint8_t data[]={1,1}; BitVector syn={2,data};
    double values[3][MAX_N]={{0.2,0.05},{0.2,0.05},{0.2,0.05}};
    QbpResult post=posterior(values,2);
    for (unsigned dai=0;dai<2;++dai) {
        QuantumLcgcdConfig cfg=config(0);
        cfg.stop_policy=dai ? QUANTUM_STOP_DAI : QUANTUM_STOP_TRIVIAL;
        QuantumDecoderWorkspace w={0};
        CHECK(quantum_decoder_init(&w,&code,&cfg)==DECODE_SUCCESS);
        QuantumDecodeResult out;
        reset_log();
        CHECK(run_post(&w,&post,&syn,&out)==DECODE_SUCCESS);
        CHECK(out.candidate_valid && out.syndrome_matched);
        verify_pauli_bound(&f,&post,3,&w,&out);
        if (dai) {
            /* Both re-encoded bits belong to qubit 0: one Eq. (25) term. */
            double r=sqrt(0.2),tau=r/(1+exp(r));
            CHECK(near(out.dai_threshold,tau));
            CHECK(out.early_stop==QUANTUM_EARLY_DAI && !out.optimality_proven);
            CHECK(out.stop_reason==DECODE_STOP_NONE && out.stats.candidates_evaluated==3);
            CHECK(out.remaining_lower_bound<out.cost && out.remaining_lower_bound+out.dai_threshold>=out.cost);
            CHECK(near(out.cost,0.2)); /* The skipped fourth pattern has cost .05. */
            double best; CHECK(oracle_best(&f,&post,3,&best)==4 && near(best,0.05) && out.cost>best);
        } else CHECK(out.optimality_proven && near(out.cost,0.05) && out.stats.candidates_evaluated==4);
        quantum_decoder_free(&w);
    }
    /* A raw binary bound would skip the final Y candidate; Pauli correction must not. */
    for (size_t p=0;p<3;++p) values[p][0]=0.055;
    values[0][1]=values[2][1]=0.4; values[1][1]=0.05;
    QuantumDecoderWorkspace w={0}; QuantumLcgcdConfig cfg=config(0);
    cfg.stop_policy=QUANTUM_STOP_TRIVIAL;
    CHECK(quantum_decoder_init(&w,&code,&cfg)==DECODE_SUCCESS);
    QuantumDecodeResult out;
    CHECK(run_post(&w,&post,&syn,&out)==DECODE_SUCCESS);
    CHECK(w.reliability[1]+w.reliability[3]>0.055);
    CHECK(out.optimality_proven && near(out.cost,0.05) && out.stats.candidates_evaluated==4);
    verify_pauli_bound(&f,&post,3,&w,&out);
    quantum_decoder_free(&w);
    /* Zero full cost proves optimal before a budget/resource check or expansion. */
    data[0]=data[1]=0;
    cfg=(QuantumLcgcdConfig){.candidate_limit=1,.slva_max_candidates=1}; /* Default trivial. */
    CHECK(quantum_decoder_init(&w,&code,&cfg)==DECODE_SUCCESS);
    CHECK(run_post(&w,&post,&syn,&out)==DECODE_SUCCESS);
    CHECK(out.cost==0 && out.optimality_proven && out.early_stop==QUANTUM_EARLY_TRIVIAL);
    CHECK(out.stats.candidates_evaluated==1 && out.stats.nodes_expanded==0);
    quantum_decoder_free(&w);
    /* Eq. (26) ties yield r=0; very large r must not overflow exp(r). */
    cfg=config(0); cfg.stop_policy=QUANTUM_STOP_DAI;
    CHECK(quantum_decoder_init(&w,&code,&cfg)==DECODE_SUCCESS);
    values[0][0]=0; values[1][0]=values[2][0]=0.2;
    CHECK(run_post(&w,&post,&syn,&out)==DECODE_SUCCESS && out.dai_threshold==0);
    for (size_t p=0;p<3;++p) values[p][0]=1e6;
    CHECK(run_post(&w,&post,&syn,&out)==DECODE_SUCCESS && out.dai_threshold==0);
    quantum_decoder_free(&w);
    quantum_code_free(&code);
    /* Mixed X/Z membership: project Ac onto qubits, with no duplicate Y term. */
    code=make_code(&toys[1]);
    CHECK(quantum_decoder_init(&w,&code,&cfg)==DECODE_SUCCESS);
    for (size_t p=0;p<3;++p) { values[p][0]=0.2; values[p][1]=0.05; }
    uint8_t redundant[]={0,0,0}; syn=(BitVector){3,redundant};
    CHECK(run_post(&w,&post,&syn,&out)==DECODE_SUCCESS);
    double r0=sqrt(0.2),r1=sqrt(0.05);
    CHECK(near(out.dai_threshold,r0/(1+exp(r0))+r1/(1+exp(r1))));
    quantum_decoder_free(&w);
    cfg=config(2); cfg.stop_policy=QUANTUM_STOP_DAI; /* Ac empty. */
    CHECK(quantum_decoder_init(&w,&code,&cfg)==DECODE_SUCCESS);
    CHECK(run_post(&w,&post,&syn,&out)==DECODE_SUCCESS && out.dai_threshold==0);
    CHECK(out.early_stop!=QUANTUM_EARLY_DAI && out.optimality_proven);
    quantum_decoder_free(&w); quantum_code_free(&code);
    CHECK(bound_stops>0 && bound_saved>0);
    fprintf(stderr,"Pauli-bound oracle: early_stops=%u saved_full_candidates=%u; DAI/non-proof counterexample passed\n",
            bound_stops,bound_saved);
}

static void test_y_metric_and_limits(void)
{
    const Fixture f={1,1,"Z"};
    QuantumCode code=make_code(&f);
    double values[3][MAX_N]={{log(10)},{log(1.25)},{log(10)}};
    QbpResult post=posterior(values,1);
    uint8_t bit=1; BitVector syn={1,&bit};
    for (size_t delta=0;delta<=1;++delta) {
        for (unsigned mode=0;mode<4;++mode) {
            QuantumDecoderWorkspace w={0};
            QuantumLcgcdConfig cfg=config(delta);
            if (mode==1) cfg.candidate_limit=1;
            if (mode==2) cfg.slva_max_candidates=1;
            if (mode==3) cfg.candidate_limit=2;
            CHECK(quantum_decoder_init(&w,&code,&cfg)==DECODE_SUCCESS);
            for (unsigned replay=0;replay<2;++replay) {
                reset_log();
                QuantumDecodeResult out;
                DecodeStatus status=run_post(&w,&post,&syn,&out);
                CHECK(out.candidate_valid && out.syndrome_matched);
                if (mode==1 || mode==2) {
                    CHECK(status==DECODE_TRUNCATED && !out.optimality_proven);
                    CHECK(out.recovery[0]==PAULI_X && near(out.cost,log(10)));
                    CHECK(out.stats.candidates_evaluated==1);
                    CHECK(out.stop_reason==(mode==1 ? DECODE_STOP_BUDGET : DECODE_STOP_RESOURCE));
                } else {
                    CHECK(status==DECODE_SUCCESS && out.optimality_proven);
                    CHECK(out.recovery[0]==PAULI_Y && near(out.cost,log(1.25)));
                    CHECK(out.stats.candidates_evaluated==2 && out.stop_reason==DECODE_STOP_EXHAUSTED);
                }
            }
            quantum_decoder_free(&w);
        }
    }
    QuantumDecoderWorkspace w={0};
    QuantumLcgcdConfig cfg=config(1); cfg.state_limit=1;
    CHECK(quantum_decoder_init(&w,&code,&cfg)==DECODE_TRUNCATED && !w.code);
    cfg=config(1); cfg.slva_max_bytes=1;
    CHECK(quantum_decoder_init(&w,&code,&cfg)==DECODE_TRUNCATED && !w.code);
    quantum_code_free(&code);
}

static void test_logical_equivalence(void)
{
    unsigned degenerate=0,logical_fail=0,syndrome_fail=0;
    for (size_t fixture=0;fixture<sizeof(toys)/sizeof(*toys);++fixture) {
        const Fixture *f=&toys[fixture];
        QuantumCode code=make_code(f);
        QuantumDecoderWorkspace w={0}; QuantumLcgcdConfig cfg=config(0);
        CHECK(quantum_decoder_init(&w,&code,&cfg)==DECODE_SUCCESS);
        uint8_t member[PATTERNS]; oracle_group(f,member);
        unsigned physical_word=(1u<<(2*f->n))-1;
        Pauli physical[MAX_N],recovery[MAX_N]; unpack(physical_word,physical,f->n);
        for (unsigned residual=0;residual<(1u<<(2*f->n));++residual) {
            unpack(physical_word^residual,recovery,f->n);
            QuantumRecoveryCheck check;
            track(1);
            QuantumStatus status=quantum_recovery_check(&w,physical,recovery,f->n,&check);
            track(0);
            CHECK(status==QUANTUM_OK);
            bool same_syn=oracle_syndrome(f,physical_word)==oracle_syndrome(f,physical_word^residual);
            CHECK(check.syndrome_matched==same_syn && check.logical_success==(member[residual]!=0));
            if (residual && check.logical_success) ++degenerate;
            if (check.syndrome_matched && !check.logical_success) ++logical_fail;
            if (!check.syndrome_matched) ++syndrome_fail;
            ++logical_cases;
        }
        quantum_decoder_free(&w); quantum_code_free(&code);
    }
    CHECK(degenerate>0 && logical_fail>0 && syndrome_fail>0);
    fprintf(stderr,"Logical oracle: cases=%u degenerate=%u syndrome-correct logical-fail=%u syndrome-fail=%u\n",
            logical_cases,degenerate,logical_fail,syndrome_fail);
}

static void test_end_to_end(void)
{
    QuantumCode code=make_code(&five);
    QuantumDecoderWorkspace failure={0},always={0};
    QuantumLcgcdConfig cfg=config(2);
    cfg.stop_policy=QUANTUM_STOP_TRIVIAL;
    CHECK(quantum_decoder_init(&failure,&code,&cfg)==DECODE_SUCCESS);
    cfg.post_policy=BP_POST_ALWAYS;
    CHECK(quantum_decoder_init(&always,&code,&cfg)==DECODE_SUCCESS);
    PauliProbabilities prob;
    double xyz[3],values[3][MAX_N];
    CHECK(quantum_depolarizing_probabilities(0.12,&prob)==QUANTUM_OK);
    CHECK(quantum_prior_llr(&prob,xyz)==QUANTUM_OK);
    for (size_t p=0;p<3;++p)
        for (size_t i=0;i<5;++i) values[p][i]=xyz[p];
    QbpPrior prior={5,values[0],values[1],values[2]};
    QbpConfig bp={.max_iterations=1,.nms_scale=0.625};
    uint8_t member[PATTERNS]; oracle_group(&five,member);
    unsigned post_frames=0,bp_frames=0,logical_errors=0,changed=0;
    uint64_t candidates=0,iterations=0;
    for (uint64_t frame=0;frame<64;++frame) {
        reset_log();
        Pauli error[5],replay_error[5]; uint8_t syn_data[4]; BitVector syn={4,syn_data};
        CHECK(quantum_depolarizing_sample(0.12,2026,frame,5,error)==QUANTUM_OK);
        CHECK(quantum_code_syndrome(&code,error,5,&syn)==QUANTUM_OK);
        unsigned target=oracle_syndrome(&five,pack(error,5));
        for (size_t r=0;r<4;++r) CHECK(syn_data[r]==((target>>r)&1u));
        QuantumDecodeResult out;
        CHECK(run_bp(&failure,&bp,&prior,&syn,&out)==DECODE_SUCCESS);
        CHECK(out.candidate_valid && out.syndrome_matched);
        CHECK(out.postprocess_invoked==!out.posterior.syndrome_matched);
        CHECK(oracle_syndrome(&five,pack(out.recovery,5))==target);
        double best;
        size_t count=oracle_best(&five,&out.posterior,target,&best);
        CHECK(count==64 && near(out.cost,best));
        if (out.postprocess_invoked) {
            ++post_frames;
            CHECK(out.stats.candidates_evaluated<=count && out.optimality_proven);
        } else {
            ++bp_frames;
            CHECK(out.stats.candidates_evaluated==0 && !out.optimality_proven);
            CHECK(memcmp(out.recovery,out.posterior.hard_estimate,5*sizeof(Pauli))==0);
        }
        candidates += out.stats.candidates_evaluated;
        iterations += out.posterior.iterations_used;
        unsigned any_changed=0;
        for (size_t i=0;i<5;++i) {
            double p[] = {0,out.posterior.post_x[i],out.posterior.post_y[i],out.posterior.post_z[i]};
            double minimum=fmin(0,fmin(p[1],fmin(p[2],p[3])));
            for (size_t a=0;a<4;++a) CHECK(near(failure.swt[i].w[a],p[a]-minimum));
            any_changed |= !near(out.posterior.post_x[i],values[0][i]);
        }
        changed += any_changed;
        QuantumRecoveryCheck logical;
        track(1);
        QuantumStatus logical_status=quantum_recovery_check(&failure,error,out.recovery,5,&logical);
        track(0);
        CHECK(logical_status==QUANTUM_OK && logical.syndrome_matched);
        CHECK(logical.logical_success==(member[pack(error,5)^pack(out.recovery,5)]!=0));
        logical_errors += !logical.logical_success;
        /* Repeat the full frame on the same workspace after logical scratch use. */
        Pauli saved[5]; memcpy(saved,out.recovery,sizeof(saved));
        double saved_post[3][5];
        memcpy(saved_post[0],out.posterior.post_x,sizeof(saved_post[0]));
        memcpy(saved_post[1],out.posterior.post_y,sizeof(saved_post[1]));
        memcpy(saved_post[2],out.posterior.post_z,sizeof(saved_post[2]));
        double saved_cost=out.cost; DecoderStats saved_stats=out.stats;
        CHECK(quantum_depolarizing_sample(0.12,2026,frame,5,replay_error)==QUANTUM_OK);
        CHECK(memcmp(error,replay_error,sizeof(error))==0);
        QuantumDecodeResult replay;
        CHECK(run_bp(&failure,&bp,&prior,&syn,&replay)==DECODE_SUCCESS);
        CHECK(memcmp(saved,replay.recovery,sizeof(saved))==0 && saved_cost==replay.cost);
        CHECK(memcmp(saved_post[0],replay.posterior.post_x,sizeof(saved_post[0]))==0);
        CHECK(memcmp(saved_post[1],replay.posterior.post_y,sizeof(saved_post[1]))==0);
        CHECK(memcmp(saved_post[2],replay.posterior.post_z,sizeof(saved_post[2]))==0);
        CHECK(memcmp(&saved_stats,&replay.stats,sizeof(saved_stats))==0);
        QuantumDecodeResult other;
        CHECK(run_bp(&always,&bp,&prior,&syn,&other)==DECODE_SUCCESS);
        CHECK(other.postprocess_invoked && other.optimality_proven && other.stats.candidates_evaluated<=count);
        CHECK(near(other.cost,best) && other.posterior.iterations_used==out.posterior.iterations_used);
        CHECK(memcmp(saved_post[0],other.posterior.post_x,sizeof(saved_post[0]))==0);
        CHECK(memcmp(saved_post[1],other.posterior.post_y,sizeof(saved_post[1]))==0);
        CHECK(memcmp(saved_post[2],other.posterior.post_z,sizeof(saved_post[2]))==0);
        if (frame==7) fprintf(stderr,"Single-frame seed=2026 frame=7: BP matched=%d post=%d candidates=%llu recovery=%u%u%u%u%u logical=%d\n",
            replay.posterior.syndrome_matched,replay.postprocess_invoked,
            (unsigned long long)replay.stats.candidates_evaluated,replay.recovery[0],replay.recovery[1],replay.recovery[2],
            replay.recovery[3],replay.recovery[4],logical.logical_success);
    }
    CHECK(post_frames>0 && bp_frames>0 && changed>0);
    for (size_t p=0;p<3;++p)
        for (size_t i=0;i<5;++i) CHECK(values[p][i]==xyz[p]);
    fprintf(stderr,"BP-LCGCD smoke: p=.12 seed=2026 frames=64 bp_direct=%u post=%u iterations=%llu candidates=%llu logical_errors=%u posterior_changed=%u\n",
        bp_frames,post_frames,(unsigned long long)iterations,(unsigned long long)candidates,logical_errors,changed);
    quantum_decoder_free(&always); quantum_decoder_free(&failure); quantum_code_free(&code);
}

int main(int argc,char **argv)
{
    if (argc>1) { debug_log=argv[1]; reset_log(); }
    test_soft_information();
    test_exhaustive_postprocessing();
    test_y_metric_and_limits();
    test_stopping_policies();
    test_logical_equivalence();
    test_end_to_end();
#ifdef QBP_TRACK_ALLOCATIONS
    CHECK(allocations==0 && frees==0 && eliminations==0);
    fprintf(stderr,"Frame/search/logical allocator calls=%u free=%u GE=%u\n",allocations,frees,eliminations);
#endif
    fprintf(stderr,"Quantum decoder passed: %u checks, %u exhaustive configurations, %u logical residuals\n",
            checks,oracle_cases,logical_cases);
    return EXIT_SUCCESS;
}
