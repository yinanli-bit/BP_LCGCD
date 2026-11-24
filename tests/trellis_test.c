#include "bplcgcd/trellis.h"

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

typedef struct {
    size_t order[7], masks[7], columns[7];
    uint8_t pa[21], pac[21], sa[3], sac[3], scratch[3];
    LocalConstraints local;
    SyndromeTrellis trellis;
} Buffers;

static void prepare(const BinarySystematic *form, const double *weights,
                    const BitVector *residual, size_t delta, Buffers *b)
{
    memset(b, 0, sizeof(*b));
    size_t m = form->rank, k = form->reduced.cols-m;
    b->local = (LocalConstraints){.basis_order=m ? b->order : NULL,
        .p_a={delta,k,delta && k ? b->pa : NULL}, .p_ac={m-delta,k,m>delta && k ? b->pac : NULL},
        .s_a={delta,delta ? b->sa : NULL}, .s_ac={m-delta,m>delta ? b->sac : NULL}};
    BitVector scratch = {form->reduced.rows, b->scratch};
    CHECK(local_constraints_prepare(form, weights, form->reduced.cols, residual, delta, &b->local, &scratch) == GF2_OK);
    b->trellis = (SyndromeTrellis){.length=delta+k,
        .column_masks=delta+k ? b->masks : NULL, .original_columns=delta+k ? b->columns : NULL};
    CHECK(syndrome_trellis_prepare(form, &b->local, weights, form->reduced.cols, 8, &b->trellis) == TRELLIS_OK);
}

static void verify(const BinarySystematic *form, const double *weights,
                   const BitVector *residual, size_t delta)
{
    Buffers b;
    prepare(form, weights, residual, delta, &b);
    SyndromeTrellis *t = &b.trellis;
    size_t m = form->rank, k = form->reduced.cols-m;
    CHECK(t->state_count == ((size_t)1 << delta) && t->flip_weights == weights);
    for (size_t stage = 0; stage < t->length; ++stage) {
        size_t original = stage < delta ? b.order[stage] : form->sys_to_orig[m+stage-delta];
        CHECK(t->original_columns[stage] == original);
        for (size_t state = 0; state < t->state_count; ++state) {
            for (uint8_t bit = 0; bit < 2; ++bit) {
                size_t next, expected = 0;
                double cost;
                CHECK(syndrome_trellis_branch(t, stage, state, bit, &next, &cost) == TRELLIS_OK);
                /* Vector parity oracle, independent of packed column masks. */
                for (size_t i = 0; i < delta; ++i) {
                    unsigned entry = stage < delta ? (i == stage) : b.pa[i*k+stage-delta];
                    expected |= (((state >> i) & 1u) ^ (bit * entry)) << i;
                }
                CHECK(next == expected && cost == bit*weights[original]);
            }
        }
    }
    unsigned accepted = 0;
    /* Test oracle only: all labeled paths, no Viterbi/SLVA production search. */
    for (unsigned path = 0; path < (1u << t->length); ++path) {
        size_t state = 0;
        double metric = 0, direct = 0;
        for (size_t stage = 0; stage < t->length; ++stage) {
            uint8_t bit = (uint8_t)((path >> stage) & 1u);
            double cost;
            CHECK(syndrome_trellis_branch(t, stage, state, bit, &state, &cost) == TRELLIS_OK);
            metric += cost;
            size_t original = stage < delta ? b.order[stage] : form->sys_to_orig[m+stage-delta];
            direct += bit*weights[original];
        }
        int valid = 1;
        for (size_t i = 0; i < delta; ++i) {
            unsigned sum = (path >> i) & 1u;
            for (size_t j = 0; j < k; ++j) sum += b.pa[i*k+j] * ((path >> (delta+j)) & 1u);
            valid &= sum%2 == b.sa[i];
        }
        CHECK((state == t->target) == valid && metric == direct);
        accepted += valid;
    }
    CHECK(accepted == (1u << k)); /* Exactly one f_A for every f_P. */
    ++graphs;
}

static void test_manual_and_errors(void)
{
    uint8_t data[] = {0,1,1,0,1,0, 0,0,1,1,0,1, 0,1,0,1,1,1};
    BinaryMatrix h = {3,6,data};
    BinarySystematic form = {0};
    CHECK(binary_matrix_systematize(&h, &form) == GF2_OK);
    double weights[] = {1000,3,6,9,8,7};
    uint8_t residual_data[] = {1,1,0};
    BitVector residual = {3,residual_data};
    for (size_t delta = 0; delta <= 2; ++delta) verify(&form, weights, &residual, delta);
    Buffers b;
    for (size_t delta = 1; delta <= 2; ++delta) {
        prepare(&form, weights, &residual, delta, &b);
        const size_t masks1[] = {1,0,1,0,1}, masks2[] = {1,2,0,3,2,3};
        const size_t columns1[] = {2,0,3,4,5}, columns2[] = {2,1,0,3,4,5};
        CHECK(!memcmp(b.masks, delta==1 ? masks1 : masks2, b.trellis.length*sizeof(size_t)));
        CHECK(!memcmp(b.columns, delta==1 ? columns1 : columns2, b.trellis.length*sizeof(size_t)));
        CHECK(b.trellis.target == 1);
        printf("Hand delta=%zu, target=1: stage original bit before after edge_cost total\n", delta);
        size_t state = 0;
        double metric = 0;
        for (size_t stage = 0; stage < b.trellis.length; ++stage) {
            uint8_t bit = (uint8_t)((10u >> b.columns[stage]) & 1u); /* f(original)=010100. */
            size_t before = state;
            double cost;
            CHECK(syndrome_trellis_branch(&b.trellis, stage, state, bit, &state, &cost) == TRELLIS_OK);
            metric += cost;
            printf("  %zu %zu %u %zu %zu %.0f %.0f\n", stage, b.columns[stage], (unsigned)bit, before, state, cost, metric);
        }
        CHECK(state == b.trellis.target && metric == (delta==1 ? 9 : 12));
    }
    size_t next;
    double cost;
    CHECK(syndrome_trellis_branch(NULL,0,0,0,&next,&cost) == TRELLIS_INVALID_ARGUMENT);
    CHECK(syndrome_trellis_branch(&b.trellis,0,0,2,&next,&cost) == TRELLIS_INVALID_ARGUMENT);
    CHECK(syndrome_trellis_branch(&b.trellis,0,0,0,NULL,&cost) == TRELLIS_INVALID_ARGUMENT);
    CHECK(syndrome_trellis_branch(&b.trellis,b.trellis.length,0,0,&next,&cost) == TRELLIS_DIMENSION_MISMATCH);
    CHECK(syndrome_trellis_branch(&b.trellis,0,4,0,&next,&cost) == TRELLIS_DIMENSION_MISMATCH);
    CHECK(syndrome_trellis_prepare(&form,&b.local,weights,6,3,&b.trellis) == TRELLIS_RESOURCE_LIMIT);
    CHECK(syndrome_trellis_prepare(&form,&b.local,weights,6,0,&b.trellis) == TRELLIS_INVALID_ARGUMENT);
    CHECK(syndrome_trellis_prepare(&form,&b.local,weights,5,4,&b.trellis) == TRELLIS_DIMENSION_MISMATCH);
    --b.trellis.length;
    CHECK(syndrome_trellis_prepare(&form,&b.local,weights,6,4,&b.trellis) == TRELLIS_DIMENSION_MISMATCH);
    ++b.trellis.length;
    const double bad[] = {-1,NAN,INFINITY};
    for (size_t i = 0; i < 3; ++i) {
        weights[2] = bad[i];
        CHECK(syndrome_trellis_prepare(&form,&b.local,weights,6,4,&b.trellis) == TRELLIS_INVALID_ARGUMENT);
    }
    weights[2] = weights[1] = DBL_MAX;
    CHECK(syndrome_trellis_prepare(&form,&b.local,weights,6,4,&b.trellis) == TRELLIS_NUMERIC_ERROR);
    binary_systematic_free(&form);
}

static void test_state_width_limit(void)
{
    size_t width = sizeof(size_t)*CHAR_BIT, n = width+1;
    BinaryMatrix h = {0};
    BinarySystematic form = {0};
    CHECK(binary_matrix_init(&h,width,n) == GF2_OK);
    for (size_t i = 0; i < width; ++i) h.data[i*n+i]=1;
    CHECK(binary_matrix_systematize(&h,&form) == GF2_OK);
    double *weights = malloc(n*sizeof(*weights));
    size_t *order = malloc(width*sizeof(*order)), *masks = malloc(n*sizeof(*masks)), *columns = malloc(n*sizeof(*columns));
    uint8_t *pa = calloc(width,1), *sa = calloc(width,1), *scratch = calloc(width,1), *residual = calloc(width,1);
    CHECK(weights && order && masks && columns && pa && sa && scratch && residual);
    for (size_t i = 0; i < n; ++i) weights[i]=1;
    LocalConstraints local = {.basis_order=order, .p_a={width,1,pa}, .p_ac={0,1,NULL}, .s_a={width,sa}};
    BitVector s = {width,residual}, work = {width,scratch};
    CHECK(local_constraints_prepare(&form,weights,n,&s,width,&local,&work) == GF2_OK);
    SyndromeTrellis t = {.length=n,.column_masks=masks,.original_columns=columns};
    CHECK(syndrome_trellis_prepare(&form,&local,weights,n,SIZE_MAX,&t) == TRELLIS_RESOURCE_LIMIT);
    uint8_t ac, sac;
    local.p_a.rows = local.s_a.length = width-1;
    local.p_ac = (BinaryMatrix){1,1,&ac};
    local.s_ac = (BitVector){1,&sac};
    CHECK(local_constraints_prepare(&form,weights,n,&s,width-1,&local,&work) == GF2_OK);
    t.length = width;
    CHECK(syndrome_trellis_prepare(&form,&local,weights,n,8,&t) == TRELLIS_RESOURCE_LIMIT);
    free(residual); free(scratch); free(sa); free(pa); free(columns); free(masks); free(order); free(weights);
    binary_systematic_free(&form);
    binary_matrix_free(&h);
}

int main(void)
{
    test_manual_and_errors();
    test_state_width_limit();
    uint32_t seed = 42;
    for (unsigned trial = 0; trial < 16; ++trial) {
        uint8_t data[21] = {0}, residual_data[3] = {0};
        BinaryMatrix h = {1+trial%3,1+trial%7,data};
        if (trial == 1) { h.rows=h.cols=2; data[0]=data[3]=1; } /* K=0, including L=0. */
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
    printf("Trellis checks passed (%u checks, %u graphs; labeled paths, Eq.9, cost and resource limits)\n",checks,graphs);
    return EXIT_SUCCESS;
}
