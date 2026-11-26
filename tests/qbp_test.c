#include "bplcgcd/qbp.h"
#include "bplcgcd/quantum_channel.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static unsigned checks;
#define CHECK(c) do { if (!(c)) { \
    fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, #c); exit(EXIT_FAILURE); \
} ++checks; } while (0)

#ifdef QBP_TRACK_ALLOCATIONS
static int track_allocations;
static unsigned hot_allocations, hot_frees;
void *__real_malloc(size_t size);
void *__real_calloc(size_t count, size_t size);
void *__real_realloc(void *ptr, size_t size);
void __real_free(void *ptr);
void *__wrap_malloc(size_t size)
{
    if (track_allocations) ++hot_allocations;
    return __real_malloc(size);
}
void *__wrap_calloc(size_t count, size_t size)
{
    if (track_allocations) ++hot_allocations;
    return __real_calloc(count, size);
}
void *__wrap_realloc(void *ptr, size_t size)
{
    if (track_allocations) ++hot_allocations;
    return __real_realloc(ptr, size);
}
void __wrap_free(void *ptr)
{
    if (track_allocations) ++hot_frees;
    __real_free(ptr);
}
#endif

static int near(double a, double b)
{
    return fabs(a-b) <= 2e-11 * fmax(1, fabs(b));
}

enum { MAX_N = 5, MAX_M = 4 };
static const unsigned anti[4][4] = {{0,0,0,0}, {0,0,1,1}, {0,1,0,1}, {0,1,1,0}};

typedef struct {
    double posterior[3][MAX_N];
    double v2c[MAX_M][MAX_N][3], c2v[MAX_M][MAX_N];
    Pauli hard[MAX_N];
    size_t iterations;
    int matched;
} DenseResult;

/* Independent deliberately dense O(degree^2) reference. Probability masses
 * implement commutation directly; no production graph, min1/min2, hard,
 * syndrome or commutation helpers. Extrinsic sums omit the destination row. */
static DenseResult dense_reference(const Pauli *s, size_t m, size_t n,
                                  double prior[3][MAX_N],
                                  const uint8_t *measured, QbpConfig cfg)
{
    DenseResult out = {0};
    double clip = cfg.message_clip ? cfg.message_clip : 50;
    for (size_t r = 0; r < m; ++r)
        for (size_t c = 0; c < n; ++c)
            for (size_t p = 0; p < 3; ++p) out.v2c[r][c][p] = prior[p][c];
    for (size_t it = 0; it < cfg.max_iterations; ++it) {
        for (size_t r = 0; r < m; ++r) {
            for (size_t c = 0; c < n; ++c) {
                if (s[r*n+c] == PAULI_I) continue;
                double magnitude = INFINITY;
                size_t neighbors = 0;
                int sign = measured[r] ? -1 : 1;
                for (size_t other = 0; other < n; ++other) {
                    Pauli label = s[r*n+other];
                    if (other == c || label == PAULI_I) continue;
                    ++neighbors;
                    double masses[] = {1, exp(-out.v2c[r][other][0]),
                        exp(-out.v2c[r][other][1]), exp(-out.v2c[r][other][2])};
                    double commute = 0, anticommute = 0;
                    for (unsigned p = 0; p < 4; ++p) {
                        if (anti[p][label]) anticommute += masses[p];
                        else commute += masses[p];
                    }
                    double lambda = log(commute / anticommute);
                    if (lambda < 0) sign = -sign;
                    if (fabs(lambda) < magnitude) magnitude = fabs(lambda);
                }
                if (!neighbors) magnitude = clip;
                out.c2v[r][c] = sign * cfg.nms_scale * magnitude;
            }
        }
        for (size_t c = 0; c < n; ++c) {
            double best = 0;
            out.hard[c] = PAULI_I;
            for (unsigned p = 1; p < 4; ++p) {
                double post = prior[p-1][c];
                for (size_t r = 0; r < m; ++r)
                    if (anti[p][s[r*n+c]]) post += out.c2v[r][c];
                out.posterior[p-1][c] = post;
                if (post < best) { best = post; out.hard[c] = (Pauli)p; }
            }
        }
        out.matched = 1;
        for (size_t r = 0; r < m; ++r) {
            unsigned syndrome = 0;
            for (size_t c = 0; c < n; ++c) syndrome ^= anti[out.hard[c]][s[r*n+c]];
            if (syndrome != measured[r]) out.matched = 0;
        }
        out.iterations = it+1;
        if (out.matched || it+1 == cfg.max_iterations) return out;
        for (size_t r = 0; r < m; ++r) {
            for (size_t c = 0; c < n; ++c) {
                if (s[r*n+c] == PAULI_I) continue;
                for (unsigned p = 1; p < 4; ++p) {
                    double extrinsic = prior[p-1][c];
                    for (size_t other = 0; other < m; ++other)
                        if (other != r && anti[p][s[other*n+c]]) extrinsic += out.c2v[other][c];
                    out.v2c[r][c][p-1] = extrinsic;
                }
            }
        }
    }
    return out;
}

static QuantumCode make_code(const Pauli *s, size_t m, size_t n)
{
    uint8_t x[MAX_M*MAX_N], z[MAX_M*MAX_N];
    CHECK(m <= MAX_M && n <= MAX_N);
    for (size_t k = 0; k < m*n; ++k) {
        x[k] = s[k] == PAULI_X || s[k] == PAULI_Y;
        z[k] = s[k] == PAULI_Z || s[k] == PAULI_Y;
    }
    BinaryMatrix hx = {m,n,x}, hz = {m,n,z};
    QuantumCode code = {0};
    CHECK(quantum_code_init(&code, &hx, &hz) == QUANTUM_OK);
    return code; /* Owning result transferred to caller. */
}

static void compare_reference(QbpWorkspace *w, const QbpResult *result,
                              const DenseResult *expected)
{
    CHECK(result->iterations_used == expected->iterations);
    CHECK(result->syndrome_matched == expected->matched);
    for (size_t c = 0; c < result->length; ++c) {
        CHECK(near(result->post_x[c], expected->posterior[0][c]));
        CHECK(near(result->post_y[c], expected->posterior[1][c]));
        CHECK(near(result->post_z[c], expected->posterior[2][c]));
        CHECK(result->hard_estimate[c] == expected->hard[c]);
    }
    for (size_t k = 0; k < w->code->graph.edge_count; ++k) {
        const QbpEdge *e = &w->code->graph.edges[k];
        CHECK(near(w->c2v[k], expected->c2v[e->row][e->col]));
        CHECK(near(w->v2c_x[k], expected->v2c[e->row][e->col][0]));
        CHECK(near(w->v2c_y[k], expected->v2c[e->row][e->col][1]));
        CHECK(near(w->v2c_z[k], expected->v2c[e->row][e->col][2]));
    }
}

static void test_commutation_and_hand_iteration(void)
{
    const double probabilities[][4] = {
        {0.7,0.1,0.1,0.1}, {0.2,0.4,0.1,0.3}, {0.1,0.1,0.7,0.1}, {0.25,0.25,0.25,0.25}
    };
    for (size_t k = 0; k < 4; ++k) {
        double xyz[3];
        for (size_t p = 0; p < 3; ++p) xyz[p] = log(probabilities[k][0]/probabilities[k][p+1]);
        for (unsigned label = 1; label < 4; ++label) {
            double commute = 0, anticommute = 0;
            for (unsigned p = 0; p < 4; ++p) {
                if (anti[p][label]) anticommute += probabilities[k][p];
                else commute += probabilities[k][p];
            }
            CHECK(near(qbp_commutation_llr(xyz, (Pauli)label), log(commute/anticommute)));
        }
    }
    const double large[] = {-1000,1000,1000};
    CHECK(near(qbp_commutation_llr(large, PAULI_X), 2000-log(2)));
    CHECK(near(qbp_commutation_llr(large, PAULI_Y), -1000));
    CHECK(near(qbp_commutation_llr(large, PAULI_Z), -1000));

    const Pauli s[] = {PAULI_X,PAULI_Z};
    QuantumCode code = make_code(s, 1, 2);
    QbpWorkspace w = {0};
    CHECK(qbp_workspace_init(&w, &code) == QUANTUM_OK);
    const double x[] = {log(4),log(8)}, y[] = {log(2),log(4)}, z[] = {log(8),log(2)};
    const QbpPrior prior = {2,x,y,z};
    uint8_t syndrome_data[] = {1};
    BitVector measured = {1,syndrome_data};
    QbpConfig cfg = {.max_iterations=10, .nms_scale=0.625};
    QbpResult result;
    CHECK(qbp_decode(&w, &cfg, &prior, &measured, &result) == QUANTUM_OK);
    CHECK(result.syndrome_matched && result.iterations_used == 1);
    CHECK(near(w.commutation[0], log(2)) && near(w.commutation[1], log(4)));
    CHECK(near(w.c2v[0], -0.625*log(4)) && near(w.c2v[1], -0.625*log(2)));
    CHECK(near(result.post_x[0], x[0]));
    CHECK(near(result.post_y[0], y[0]-0.625*log(4)));
    CHECK(near(result.post_z[0], z[0]-0.625*log(4)));
    CHECK(near(result.post_x[1], x[1]-0.625*log(2)));
    CHECK(near(result.post_y[1], y[1]-0.625*log(2)));
    CHECK(near(result.post_z[1], z[1]));
    CHECK(result.hard_estimate[0] == PAULI_Y && result.hard_estimate[1] == PAULI_I);
    CHECK(x[0] == log(4) && y[0] == log(2) && z[0] == log(8)); /* Prior untouched. */
    printf("Hand iteration: syndrome=1, lambda=(ln2,ln4), C2V=(-.625ln4,-.625ln2), hard=YI, iterations=1\n");
    syndrome_data[0] = 0;
    CHECK(qbp_decode(&w, &cfg, &prior, &measured, &result) == QUANTUM_OK);
    CHECK(result.syndrome_matched && result.iterations_used == 1);
    CHECK(w.c2v[0] > 0 && w.c2v[1] > 0);
    CHECK(result.hard_estimate[0] == PAULI_I && result.hard_estimate[1] == PAULI_I);
    qbp_workspace_free(&w);
    quantum_code_free(&code);
}

static void test_dense_messages(void)
{
    const Pauli s[] = {
        PAULI_X,PAULI_Y,PAULI_Z,PAULI_I,PAULI_I,
        PAULI_Z,PAULI_X,PAULI_I,PAULI_Y,PAULI_I,
        PAULI_Y,PAULI_Z,PAULI_Z,PAULI_Y,PAULI_I,
        PAULI_I,PAULI_I,PAULI_I,PAULI_I,PAULI_I
    };
    QuantumCode code = make_code(s, 4, 5);
    QbpWorkspace w = {0}, independent = {0};
    CHECK(qbp_workspace_init(&w, &code) == QUANTUM_OK);
    CHECK(qbp_workspace_init(&independent, &code) == QUANTUM_OK);
    CHECK(w.messages != independent.messages);
    unsigned multi_iteration = 0, negative_lambda = 0;
    for (unsigned frame = 0; frame < 24; ++frame) {
        double values[3][MAX_N];
        for (size_t p = 0; p < 3; ++p)
            for (size_t c = 0; c < 5; ++c)
                values[p][c] = frame == 0 ? 0 : ((int)((frame*7+c*3+p*11)%19)-9)*0.25;
        uint8_t syn[] = {(uint8_t)(frame&1), (uint8_t)((frame>>1)&1), (uint8_t)((frame>>2)&1), 0};
        BitVector measured = {4,syn};
        QbpPrior prior = {5,values[0],values[1],values[2]};
        for (size_t iterations = 1; iterations <= 7; iterations += 3) {
            QbpConfig cfg = {.max_iterations=iterations, .nms_scale=0.625};
            QbpResult result;
            CHECK(qbp_decode(&w, &cfg, &prior, &measured, &result) == QUANTUM_OK);
            DenseResult expected = dense_reference(s, 4, 5, values, syn, cfg);
            compare_reference(&w, &result, &expected);
            multi_iteration += result.iterations_used > 1;
            for (size_t k = 0; k < code.graph.edge_count; ++k) negative_lambda += w.commutation[k] < 0;
        }
    }
    CHECK(multi_iteration > 0 && negative_lambda > 0);
    /* Distinct min1/min2, repeated minima, negative and exact zero lambda. */
    const Pauli row[] = {PAULI_X,PAULI_Y,PAULI_Z,PAULI_X};
    QuantumCode row_code = make_code(row, 1, 4);
    QbpWorkspace row_w = {0};
    CHECK(qbp_workspace_init(&row_w, &row_code) == QUANTUM_OK);
    for (unsigned scenario = 0; scenario < 4; ++scenario) {
        double values[3][MAX_N] = {{0}};
        for (size_t c = 0; c < 4; ++c) {
            for (size_t p = 0; p < 3; ++p)
                values[p][c] = scenario == 0 ? 1 : ((double)c-(double)p)*0.75;
        }
        if (scenario == 2) for (size_t p = 0; p < 3; ++p) values[p][0] = 0;
        if (scenario == 3) for (size_t p = 0; p < 3; ++p) values[p][0] = values[p][1] = 0;
        for (unsigned bit = 0; bit < 2; ++bit) {
            uint8_t syn[] = {(uint8_t)bit};
            BitVector measured = {1,syn};
            QbpPrior prior = {4,values[0],values[1],values[2]};
            QbpConfig cfg = {.max_iterations=1, .nms_scale=1};
            QbpResult result;
            CHECK(qbp_decode(&row_w, &cfg, &prior, &measured, &result) == QUANTUM_OK);
            DenseResult expected = dense_reference(row, 1, 4, values, syn, cfg);
            compare_reference(&row_w, &result, &expected);
        }
    }
    /* Ordinary checks remain unsaturated above 50, even with an explicit clip. */
    double strong[3][MAX_N] = {{0}};
    for (size_t p = 0; p < 3; ++p)
        for (size_t c = 0; c < 4; ++c) strong[p][c] = 80 + 20*c;
    const double clips[] = {0,1};
    for (size_t k = 0; k < 2; ++k) {
        for (unsigned bit = 0; bit < 2; ++bit) {
            uint8_t syn[] = {(uint8_t)bit};
            BitVector measured = {1,syn};
            QbpPrior prior = {4,strong[0],strong[1],strong[2]};
            QbpConfig cfg = {.max_iterations=1, .nms_scale=0.625, .message_clip=clips[k]};
            QbpResult result;
            CHECK(qbp_decode(&row_w, &cfg, &prior, &measured, &result) == QUANTUM_OK);
            DenseResult expected = dense_reference(row, 1, 4, strong, syn, cfg);
            compare_reference(&row_w, &result, &expected);
            double sign = bit ? -1 : 1;
            CHECK(near(row_w.c2v[0], sign*0.625*(100-log(2))));
            for (size_t c = 1; c < 4; ++c)
                CHECK(near(row_w.c2v[c], sign*0.625*(80-log(2))));
        }
    }
    qbp_workspace_free(&row_w);
    quantum_code_free(&row_code);
    qbp_workspace_free(&independent);
    qbp_workspace_free(&w);
    quantum_code_free(&code);
}

static void test_empty_singleton_and_ties(void)
{
    const Pauli empty[] = {PAULI_I};
    QuantumCode code = make_code(empty, 1, 1);
    QbpWorkspace w = {0};
    CHECK(qbp_workspace_init(&w, &code) == QUANTUM_OK);
    const double cases[5][3] = {{0,-0.0,0}, {-1,-1,0}, {2,-2,-1}, {0,0,-3}, {1,2,3}};
    const Pauli expected[] = {PAULI_I,PAULI_X,PAULI_Y,PAULI_Z,PAULI_I};
    uint8_t syn_data[] = {0};
    BitVector measured = {1,syn_data};
    QbpConfig cfg = {.max_iterations=3, .nms_scale=0.625};
    QbpResult result;
    for (size_t k = 0; k < 5; ++k) {
        QbpPrior prior = {1,&cases[k][0],&cases[k][1],&cases[k][2]};
        CHECK(qbp_decode(&w, &cfg, &prior, &measured, &result) == QUANTUM_OK);
        CHECK(result.hard_estimate[0] == expected[k] && result.syndrome_matched);
        CHECK(result.iterations_used == 1 && result.post_x[0] == cases[k][0]);
    }
    syn_data[0] = 1; /* An impossible empty-row syndrome never converges. */
    QbpPrior prior = {1,&cases[0][0],&cases[0][1],&cases[0][2]};
    CHECK(qbp_decode(&w, &cfg, &prior, &measured, &result) == QUANTUM_OK);
    CHECK(!result.syndrome_matched && result.iterations_used == 3);
    cfg.nms_scale = 0;
    CHECK(qbp_decode(&w, &cfg, &prior, &measured, &result) == QUANTUM_INVALID_ARGUMENT);
    CHECK(!result.hard_estimate && !result.post_x && !result.iterations_used);
    qbp_workspace_free(&w);
    quantum_code_free(&code);

    const Pauli singleton[] = {PAULI_Y,PAULI_I};
    code = make_code(singleton, 1, 2);
    CHECK(qbp_workspace_init(&w, &code) == QUANTUM_OK);
    double values[] = {1,1};
    prior = (QbpPrior){2,values,values,values};
    cfg = (QbpConfig){.max_iterations=3, .nms_scale=0.625, .message_clip=40};
    CHECK(qbp_decode(&w, &cfg, &prior, &measured, &result) == QUANTUM_OK);
    CHECK(result.syndrome_matched && result.iterations_used == 1);
    CHECK(w.c2v[0] == -25 && result.hard_estimate[0] == PAULI_X);
    CHECK(result.post_x[0] == -24 && result.post_y[0] == 1 && result.post_z[0] == -24);
    CHECK(result.post_x[1] == 1 && result.hard_estimate[1] == PAULI_I);
    cfg.message_clip = 0; /* Default 50 applies only to this singleton. */
    for (unsigned bit = 0; bit < 2; ++bit) {
        syn_data[0] = (uint8_t)bit;
        CHECK(qbp_decode(&w, &cfg, &prior, &measured, &result) == QUANTUM_OK);
        CHECK(w.c2v[0] == (bit ? -31.25 : 31.25));
        CHECK(result.syndrome_matched && result.iterations_used == 1);
        CHECK(result.hard_estimate[0] == (bit ? PAULI_X : PAULI_I));
    }
    values[0] = DBL_MAX;
    syn_data[0] = 0;
    cfg.nms_scale = 1;
    cfg.message_clip = DBL_MAX;
    CHECK(qbp_decode(&w, &cfg, &prior, &measured, &result) == QUANTUM_NUMERIC_ERROR);
    CHECK(!result.hard_estimate && !result.post_x);
    /* An error must not poison a later independent frame. */
    values[0] = 1;
    cfg.message_clip = 40;
    CHECK(qbp_decode(&w, &cfg, &prior, &measured, &result) == QUANTUM_OK);
    CHECK(result.syndrome_matched && result.hard_estimate[0] == PAULI_I);
    qbp_workspace_free(&w);
    quantum_code_free(&code);
}

static void test_depolarizing_smoke_replay(void)
{
    /* Generic five-qubit stabilizer example from 2025 Eq. (1), no code family API. */
    const Pauli s[] = {
        PAULI_X,PAULI_Z,PAULI_Z,PAULI_X,PAULI_I,
        PAULI_Z,PAULI_Z,PAULI_X,PAULI_I,PAULI_X,
        PAULI_Z,PAULI_X,PAULI_I,PAULI_X,PAULI_Z,
        PAULI_X,PAULI_I,PAULI_X,PAULI_Z,PAULI_Z
    };
    QuantumCode code = make_code(s, 4, 5);
    QbpWorkspace w = {0}, replay_w = {0};
    CHECK(qbp_workspace_init(&w, &code) == QUANTUM_OK);
    CHECK(qbp_workspace_init(&replay_w, &code) == QUANTUM_OK);
    PauliProbabilities prob;
    double xyz[3], values[3][MAX_N];
    CHECK(quantum_depolarizing_probabilities(0.12, &prob) == QUANTUM_OK);
    CHECK(quantum_prior_llr(&prob, xyz) == QUANTUM_OK);
    for (size_t p = 0; p < 3; ++p)
        for (size_t c = 0; c < 5; ++c) values[p][c] = xyz[p];
    QbpPrior prior = {5,values[0],values[1],values[2]};
    QbpConfig cfg = {.max_iterations=12, .nms_scale=0.625};
    unsigned converged = 0, nonzero = 0, total_iterations = 0;
    double *original_messages = w.messages;
    QbpEdge *original_edges = code.graph.edges;
    for (uint64_t frame = 0; frame < 256; ++frame) {
#ifdef QBP_TRACK_ALLOCATIONS
        track_allocations = 1;
#endif
        Pauli error[5], replay_error[5];
        uint8_t syn_data[4], recovered_syn_data[4];
        BitVector measured = {4,syn_data}, recovered = {4,recovered_syn_data};
        CHECK(quantum_depolarizing_sample(0.12, 2026, frame, 5, error) == QUANTUM_OK);
        CHECK(quantum_code_syndrome(&code, error, 5, &measured) == QUANTUM_OK);
        unsigned any = 0;
        for (size_t r = 0; r < 4; ++r) {
            unsigned dense_syn = 0;
            for (size_t c = 0; c < 5; ++c) dense_syn ^= anti[error[c]][s[r*5+c]];
            CHECK(syn_data[r] == dense_syn);
            any |= syn_data[r];
        }
        nonzero += any != 0;
        QbpResult result, replay;
        CHECK(qbp_decode(&w, &cfg, &prior, &measured, &result) == QUANTUM_OK);
        DenseResult expected = dense_reference(s, 4, 5, values, syn_data, cfg);
        compare_reference(&w, &result, &expected);
        CHECK(quantum_code_syndrome(&code, result.hard_estimate, 5, &recovered) == QUANTUM_OK);
        CHECK(result.syndrome_matched == (memcmp(syn_data, recovered_syn_data, 4) == 0));
        converged += result.syndrome_matched;
        total_iterations += (unsigned)result.iterations_used;
        /* Replay frames after intervening calls, on shared topology. */
        CHECK(quantum_depolarizing_sample(0.12, 2026, frame, 5, replay_error) == QUANTUM_OK);
        CHECK(memcmp(error, replay_error, sizeof(error)) == 0);
        CHECK(qbp_decode(&replay_w, &cfg, &prior, &measured, &replay) == QUANTUM_OK);
        CHECK(memcmp(result.hard_estimate, replay.hard_estimate, 5*sizeof(Pauli)) == 0);
        CHECK(memcmp(result.post_x, replay.post_x, 5*sizeof(double)) == 0);
        CHECK(memcmp(result.post_y, replay.post_y, 5*sizeof(double)) == 0);
        CHECK(memcmp(result.post_z, replay.post_z, 5*sizeof(double)) == 0);
        CHECK(result.iterations_used == replay.iterations_used && result.syndrome_matched == replay.syndrome_matched);
        CHECK(w.messages == original_messages && code.graph.edges == original_edges);
#ifdef QBP_TRACK_ALLOCATIONS
        track_allocations = 0;
#endif
        if (frame == 7) {
            printf("Single depolarizing frame seed=2026 index=7 p=0.12: error=");
            for (size_t c = 0; c < 5; ++c) putchar("IXYZ"[error[c]]);
            printf(" syndrome=");
            for (size_t r = 0; r < 4; ++r) printf("%u", syn_data[r]);
            printf(" hard=");
            for (size_t c = 0; c < 5; ++c) putchar("IXYZ"[result.hard_estimate[c]]);
            printf(" iterations=%zu matched=%d\n", result.iterations_used, result.syndrome_matched);
        }
    }
    CHECK(nonzero > 0 && converged > 0 && converged < 256);
#ifdef QBP_TRACK_ALLOCATIONS
    CHECK(hot_allocations == 0 && hot_frees == 0);
    printf("Frame-loop allocator calls: malloc/calloc/realloc=%u free=%u\n", hot_allocations, hot_frees);
#endif
    for (size_t p = 0; p < 3; ++p)
        for (size_t c = 0; c < 5; ++c) CHECK(values[p][c] == xyz[p]);
    printf("QNMS smoke p=0.12 seed=2026: frames=256 nonzero_syndrome=%u matched=%u iterations_total=%u (not LER)\n",
           nonzero, converged, total_iterations);
    qbp_workspace_free(&replay_w);
    qbp_workspace_free(&w);
    quantum_code_free(&code);
}

int main(void)
{
    test_commutation_and_hand_iteration();
    test_dense_messages();
    test_empty_singleton_and_ties();
    test_depolarizing_smoke_replay();
    printf("QNMS passed (%u checks; dense message oracle, hand iteration, replay)\n", checks);
    return EXIT_SUCCESS;
}
