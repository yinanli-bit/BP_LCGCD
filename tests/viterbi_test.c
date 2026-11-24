#include "bplcgcd/viterbi.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static unsigned checks, graphs;
#define CHECK(condition) do { \
    if (!(condition)) { \
        fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, #condition); \
        exit(EXIT_FAILURE); \
    } \
    ++checks; \
} while (0)

/* Independent vector parity for a labeled prefix; no production masks/edges. */
static unsigned prefix_state(const LocalConstraints *local, unsigned prefix, size_t layer)
{
    size_t delta = local->p_a.rows, k = local->p_a.cols;
    unsigned state = 0;
    for (size_t i = 0; i < delta; ++i) {
        unsigned sum = i < layer ? ((prefix >> i) & 1u) : 0;
        for (size_t j = 0; j < k && delta+j < layer; ++j) sum += local->p_a.data[i*k+j] * ((prefix >> (delta+j)) & 1u);
        state |= (sum%2) << i;
    }
    return state;
}

static void verify(const BinarySystematic *form, const double *weights, const BitVector *residual, size_t delta)
{
    size_t m = form->rank, k = form->reduced.cols-m;
    size_t order[7], masks[7], columns[7];
    uint8_t pa[21], pac[21], sa[3], sac[3], scratch_data[3], bits[7], repeat_bits[7];
    LocalConstraints local = {.basis_order=m ? order : NULL,
        .p_a={delta,k,delta && k ? pa : NULL}, .p_ac={m-delta,k,m>delta && k ? pac : NULL},
        .s_a={delta,delta ? sa : NULL}, .s_ac={m-delta,m>delta ? sac : NULL}};
    BitVector scratch = {form->reduced.rows,scratch_data};
    CHECK(local_constraints_prepare(form,weights,form->reduced.cols,residual,delta,&local,&scratch) == GF2_OK);
    SyndromeTrellis t = {.length=delta+k,.column_masks=masks,.original_columns=columns};
    CHECK(syndrome_trellis_prepare(form,&local,weights,form->reduced.cols,8,&t) == TRELLIS_OK);
    ViterbiCell cells[64];
    BitVector path = {t.length,t.length ? bits : NULL};
    double cost;
    CHECK(viterbi_decode(&t,cells,64,&path,&cost) == TRELLIS_OK);
    unsigned expected_path = 0;
    for (size_t layer = 0; layer <= t.length; ++layer) {
        double optimum[8];
        unsigned chosen[8] = {0};
        for (size_t state = 0; state < t.state_count; ++state) optimum[state]=INFINITY;
        /* Exhaustive prefix oracle checks EVERY cell, including unreachable
         * states and exact ties (integer/binary-fraction weights). */
        for (unsigned prefix = 0; prefix < (1u << layer); ++prefix) {
            unsigned state = prefix_state(&local,prefix,layer);
            double metric = 0;
            for (size_t stage = 0; stage < layer; ++stage) {
                size_t original = stage < delta ? order[stage] : form->sys_to_orig[m+stage-delta];
                metric += ((prefix >> stage) & 1u)*weights[original];
            }
            if (metric < optimum[state]) { optimum[state]=metric; chosen[state]=prefix; }
        }
        for (size_t state = 0; state < t.state_count; ++state) {
            const ViterbiCell *cell = &cells[layer*t.state_count+state];
            CHECK(cell->metric == optimum[state]);
            if (!layer || !isfinite(optimum[state])) CHECK(cell->previous_state == SIZE_MAX);
            else {
                CHECK(cell->bit == ((chosen[state] >> (layer-1)) & 1u));
                CHECK(cell->previous_state == prefix_state(&local,chosen[state],layer-1));
            }
        }
        if (layer == t.length) {
            CHECK(cost == optimum[t.target]);
            expected_path=chosen[t.target];
        }
    }
    for (size_t stage = 0; stage < t.length; ++stage) CHECK(bits[stage] == ((expected_path >> stage) & 1u));
    BitVector repeat = {t.length,t.length ? repeat_bits : NULL};
    double repeat_cost;
    CHECK(viterbi_decode(&t,cells,64,&repeat,&repeat_cost) == TRELLIS_OK && repeat_cost == cost);
    CHECK(!memcmp(bits,repeat_bits,t.length));
    ++graphs;
}

static void test_manual(void)
{
    uint8_t data[] = {0,1,1,0,1,0, 0,0,1,1,0,1, 0,1,0,1,1,1};
    BinaryMatrix h = {3,6,data};
    BinarySystematic form = {0};
    CHECK(binary_matrix_systematize(&h,&form) == GF2_OK);
    double weights[] = {1000,3,6,9,8,7};
    uint8_t residual_data[] = {1,1,0};
    BitVector residual = {3,residual_data};
    for (size_t delta = 0; delta <= form.rank; ++delta) verify(&form,weights,&residual,delta);
    size_t masks[] = {1,2,0,3,2,3}, columns[] = {2,1,0,3,4,5};
    SyndromeTrellis t = {.length=6,.delta=2,.state_count=4,.target=1,
        .column_masks=masks,.original_columns=columns,.flip_weights=weights};
    ViterbiCell cells[28];
    uint8_t bits[6];
    BitVector path = {6,bits};
    double cost;
    CHECK(viterbi_decode(&t,cells,28,&path,&cost) == TRELLIS_OK && cost == 6);
    CHECK(!memcmp(bits,(uint8_t[]){1,0,0,0,0,0},6));
    const double expected[7][4] = {{0,INFINITY,INFINITY,INFINITY},{0,6,INFINITY,INFINITY},
        {0,6,3,9},{0,6,3,9},{0,6,3,9},{0,6,3,9},{0,6,3,7}};
    puts("Hand delta=2: layer metrics[q=0,1,2,3]");
    for (size_t layer = 0; layer <= 6; ++layer) {
        printf("  %zu:",layer);
        for (size_t state = 0; state < 4; ++state) {
            CHECK(cells[layer*4+state].metric == expected[layer][state]);
            printf(" %.0f",cells[layer*4+state].metric);
        }
        putchar('\n');
    }
    puts("Best path(stage order)=100000, f_I[A]=10, f_P=0000, local cost=6");
    binary_systematic_free(&form);
}

static void test_limits_and_no_path(void)
{
    size_t masks[] = {0,0}, columns[] = {0,1};
    double weights[] = {1,2}, cost = 99;
    SyndromeTrellis t = {.length=1,.delta=1,.state_count=2,.target=1,
        .column_masks=masks,.original_columns=columns,.flip_weights=weights};
    ViterbiCell cells[6];
    uint8_t bits[2] = {9,9};
    BitVector path = {1,bits};
    /* Connectivity probe: zero column cannot reach nonzero syndrome. */
    CHECK(viterbi_decode(&t,cells,6,&path,&cost) == TRELLIS_NO_PATH);
    CHECK(isinf(cells[3].metric) && cells[3].previous_state == SIZE_MAX);
    cells[0].metric=77;
    CHECK(viterbi_decode(&t,cells,3,&path,&cost) == TRELLIS_RESOURCE_LIMIT);
    CHECK(cells[0].metric == 77 && bits[0] == 9 && cost == 99);
    CHECK(viterbi_decode(NULL,cells,6,&path,&cost) == TRELLIS_INVALID_ARGUMENT);
    CHECK(viterbi_decode(&t,NULL,6,&path,&cost) == TRELLIS_INVALID_ARGUMENT);
    CHECK(viterbi_decode(&t,cells,6,&path,NULL) == TRELLIS_INVALID_ARGUMENT);
    path.length=2;
    CHECK(viterbi_decode(&t,cells,6,&path,&cost) == TRELLIS_DIMENSION_MISMATCH);
    t.length=path.length=SIZE_MAX;
    CHECK(viterbi_decode(&t,cells,SIZE_MAX,&path,&cost) == TRELLIS_RESOURCE_LIMIT);
    t.length=path.length=2;
    t.state_count=(size_t)1 << (sizeof(size_t)*CHAR_BIT-1);
    CHECK(viterbi_decode(&t,cells,SIZE_MAX,&path,&cost) == TRELLIS_RESOURCE_LIMIT);
    t.length=path.length=0;
    CHECK(viterbi_decode(&t,cells,SIZE_MAX,&path,&cost) == TRELLIS_RESOURCE_LIMIT); /* Byte size. */
    t.length=path.length=2;
    t.state_count=2; t.target=0; masks[0]=masks[1]=1;
    weights[0]=weights[1]=DBL_MAX; /* Finite edges, overflowing path sum. */
    CHECK(viterbi_decode(&t,cells,6,&path,&cost) == TRELLIS_NUMERIC_ERROR);
}

static void test_long_path(void)
{
    size_t masks[70] = {0}, columns[70];
    double weights[70];
    for (size_t j = 0; j < 70; ++j) { columns[j]=j; weights[j]=2; }
    masks[0]=masks[69]=1; weights[0]=9; weights[69]=0.25;
    SyndromeTrellis t = {.length=70,.delta=1,.state_count=2,.target=1,
        .column_masks=masks,.original_columns=columns,.flip_weights=weights};
    ViterbiCell cells[142];
    uint8_t bits[70];
    BitVector path = {70,bits};
    double cost;
    CHECK(viterbi_decode(&t,cells,142,&path,&cost) == TRELLIS_OK && cost == 0.25);
    for (size_t j = 0; j < 70; ++j) CHECK(bits[j] == (j==69));
}

int main(void)
{
    test_manual();
    test_limits_and_no_path();
    test_long_path();
    uint32_t seed = 2026;
    for (unsigned trial = 0; trial < 24; ++trial) {
        uint8_t data[21] = {0}, residual_data[3] = {0};
        BinaryMatrix h = {1+trial%3,1+trial%7,data};
        if (trial == 1) { h.rows=h.cols=2; data[0]=data[3]=1; }
        else if (trial%8) {
            for (size_t j = 0; j < h.rows*h.cols; ++j) { seed=1664525u*seed+1013904223u; data[j]=(uint8_t)(seed>>31); }
        }
        for (size_t i = 0; i < h.rows; ++i) {
            for (size_t j = 0; j < h.cols; ++j) residual_data[i] ^= data[i*h.cols+j] & (uint8_t)((37u>>j)&1u);
        }
        BinarySystematic form = {0};
        CHECK(binary_matrix_systematize(&h,&form) == GF2_OK);
        double weights[7];
        for (size_t j = 0; j < h.cols; ++j) weights[j] = trial%3 ? (double)((j+trial)%5)/4 : 0;
        BitVector residual = {h.rows,residual_data};
        for (size_t delta = 0; delta <= form.rank; ++delta) verify(&form,weights,&residual,delta);
        binary_systematic_free(&form);
    }
    printf("Viterbi checks passed (%u checks, %u oracle graphs; all prefix metrics, traceback, ties and resource limits)\n",checks,graphs);
    return EXIT_SUCCESS;
}
