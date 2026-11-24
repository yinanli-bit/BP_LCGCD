#include "bplcgcd/local_constraints.h"
#include "bplcgcd/decoder.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

static unsigned checks, partitions;
#define CHECK(condition) do { \
    if (!(condition)) { \
        fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, #condition); \
        exit(EXIT_FAILURE); \
    } \
    ++checks; \
} while (0)

/* Fixed test storage demonstrates caller-owned buffers, including zero shapes. */
typedef struct {
    size_t order[7];
    uint8_t a[28], ac[28], sa[4], sac[4], scratch[4];
    LocalConstraints local;
    BitVector transformed;
} Buffers;

static void buffers_init(Buffers *b, const BinarySystematic *form, size_t delta)
{
    size_t m = form->rank, k = form->reduced.cols - m;
    memset(b, 0, sizeof(*b));
    b->local = (LocalConstraints){.basis_order=m ? b->order : NULL,
        .p_a={delta, k, delta && k ? b->a : NULL},
        .p_ac={m-delta, k, m>delta && k ? b->ac : NULL},
        .s_a={delta, delta ? b->sa : NULL}, .s_ac={m-delta, m>delta ? b->sac : NULL}};
    b->transformed = (BitVector){form->reduced.rows, b->scratch};
}

static unsigned syndrome(const BinaryMatrix *h, unsigned word)
{
    unsigned result = 0;
    for (size_t i = 0; i < h->rows; ++i) {
        unsigned sum = 0;
        for (size_t j = 0; j < h->cols; ++j) sum += h->data[i * h->cols + j] * ((word >> j) & 1u);
        result |= (sum % 2) << i;
    }
    return result;
}

static void verify(const BinaryMatrix *h, const BinarySystematic *form,
                   const double *weights, const BitVector *residual, size_t delta)
{
    Buffers b;
    buffers_init(&b, form, delta);
    CHECK(local_constraints_prepare(form, weights, h->cols, residual, delta, &b.local, &b.transformed) == GF2_OK);
    size_t m = form->rank, k = h->cols - m;
    /* Independent ranking oracle counts how many BASIS columns precede each.
     * No production sorting, and free-column weights do not participate. */
    for (size_t original = 0; original < h->cols; ++original) {
        if (form->orig_to_sys[original] >= m) continue;
        size_t ahead = 0;
        for (size_t other = 0; other < h->cols; ++other) {
            if (form->orig_to_sys[other] < m && (weights[other] > weights[original] ||
                (weights[other] == weights[original] && other < original))) ++ahead;
        }
        CHECK(b.order[ahead] == original);
    }
    for (size_t part = 0; part < 2; ++part) {
        const BinaryMatrix *p = part ? &b.local.p_ac : &b.local.p_a;
        const BitVector *s = part ? &b.local.s_ac : &b.local.s_a;
        size_t start = part ? delta : 0;
        for (size_t i = 0; i < s->length; ++i) {
            size_t basis = form->orig_to_sys[b.order[start + i]];
            unsigned sum = 0;
            for (size_t r = 0; r < h->rows; ++r) sum += form->row_transform.data[basis * h->rows + r] * residual->data[r];
            CHECK(s->data[i] == sum % 2);
            for (size_t j = 0; j < k; ++j) CHECK(p->data[i*k+j] == form->reduced.data[basis*h->cols+m+j]);
        }
    }
    unsigned target = 0;
    for (size_t i = 0; i < h->rows; ++i) target |= (unsigned)residual->data[i] << i;
    /* All original-coordinate flips: Eq.9 AND Eq.10 iff H*f=s_res. */
    for (unsigned f = 0; f < (1u << h->cols); ++f) {
        int split_valid = 1;
        for (size_t part = 0; part < 2; ++part) {
            const BinaryMatrix *p = part ? &b.local.p_ac : &b.local.p_a;
            const BitVector *s = part ? &b.local.s_ac : &b.local.s_a;
            size_t start = part ? delta : 0;
            for (size_t i = 0; i < s->length; ++i) {
                unsigned sum = (f >> b.order[start + i]) & 1u;
                for (size_t j = 0; j < k; ++j) sum += p->data[i*k+j] * ((f >> form->sys_to_orig[m+j]) & 1u);
                split_valid &= sum % 2 == s->data[i];
            }
        }
        CHECK(split_valid == (syndrome(h, f) == target));
    }
    ++partitions;
}

static void test_manual_and_modes(void)
{
    uint8_t data[] = {0,1,1,0,1,0, 0,0,1,1,0,1, 0,1,0,1,1,1};
    BinaryMatrix h = {3, 6, data};
    BinarySystematic form = {0};
    CHECK(binary_matrix_systematize(&h, &form) == GF2_OK && form.rank == 2);
    CHECK(form.sys_to_orig[0] == 1 && form.sys_to_orig[1] == 2);
    uint8_t target_data[] = {1,1,0};
    BitVector target = {3, target_data};
    double weights[] = {1000,3,6,900,800,700};
    Buffers b;
    buffers_init(&b, &form, 1);
    CHECK(local_constraints_prepare(&form, weights, 6, &target, 1, &b.local, &b.transformed) == GF2_OK);
    CHECK(b.order[0] == 2 && b.order[1] == 1); /* Larger FREE weights ignored. */
    const uint8_t pa[] = {0,1,0,1}, pac[] = {0,1,1,1};
    CHECK(!memcmp(b.a, pa, 4) && !memcmp(b.ac, pac, 4));
    CHECK(b.sa[0] == 1 && b.sac[0] == 0);
    for (size_t delta = 0; delta <= form.rank; ++delta) verify(&h, &form, weights, &target, delta);
    weights[1] = weights[2];
    CHECK(local_constraints_prepare(&form, weights, 6, &target, 1, &b.local, &b.transformed) == GF2_OK);
    CHECK(b.order[0] == 1 && b.order[1] == 2); /* Ties: original index. */
    weights[1] = 3;
    double llr[] = {1000,-3,6,900,800,-700}, normalized[6];
    uint8_t hint_data[] = {0,1,0,0,0,1}, reference_data[6], residual_data[3];
    BitVector hint = {6, hint_data}, reference = {6, reference_data}, residual = {3, residual_data};
    for (unsigned mode = 0; mode < 3; ++mode) {
        DecodeInput input = {.mode=(DecodeMode)mode, .h=&h, .length=6};
        if (mode == DECODE_SYNDROME) { input.weights=weights; input.z_hint=&hint; }
        else input.llr=llr;
        if (mode != DECODE_CLASSICAL) input.s_target=&target;
        CHECK(decode_prepare(&input, &reference, &residual, normalized) == DECODE_SUCCESS);
        for (size_t delta = 0; delta <= form.rank; ++delta) verify(&h, &form, normalized, &residual, delta);
    }
    CHECK(!memcmp(data, (uint8_t[]){0,1,1,0,1,0, 0,0,1,1,0,1, 0,1,0,1,1,1}, sizeof(data)));
    CHECK(!memcmp(weights, (double[]){1000,3,6,900,800,700}, sizeof(weights)));
    CHECK(!memcmp(target_data, (uint8_t[]){1,1,0}, sizeof(target_data)));
    CHECK(local_constraints_prepare(NULL, weights, 6, &target, 1, &b.local, &b.transformed) == GF2_INVALID_ARGUMENT);
    CHECK(local_constraints_prepare(&form, weights, 5, &target, 1, &b.local, &b.transformed) == GF2_DIMENSION_MISMATCH);
    CHECK(local_constraints_prepare(&form, weights, 6, &target, 3, &b.local, &b.transformed) == GF2_DIMENSION_MISMATCH);
    b.local.p_a.cols = 3;
    CHECK(local_constraints_prepare(&form, weights, 6, &target, 1, &b.local, &b.transformed) == GF2_DIMENSION_MISMATCH);
    buffers_init(&b, &form, 1);
    b.local.p_a.data = NULL;
    CHECK(local_constraints_prepare(&form, weights, 6, &target, 1, &b.local, &b.transformed) == GF2_INVALID_ARGUMENT);
    buffers_init(&b, &form, 1);
    const double bad[] = {-1, NAN, INFINITY};
    for (size_t i = 0; i < 3; ++i) {
        weights[0] = bad[i];
        CHECK(local_constraints_prepare(&form, weights, 6, &target, 1, &b.local, &b.transformed) == GF2_INVALID_ARGUMENT);
    }
    weights[0] = 1000;
    target_data[2] ^= 1;
    CHECK(local_constraints_prepare(&form, weights, 6, &target, 1, &b.local, &b.transformed) == GF2_INCONSISTENT_SYNDROME);
    binary_systematic_free(&form);
    puts("Hand example: A(original)={2}, Ac={1}; P_A=0101, s_A=1; P_Ac=0111, s_Ac=0");
}

int main(void)
{
    test_manual_and_modes();
    uint32_t seed = 12345;
    for (unsigned trial = 0; trial < 24; ++trial) {
        uint8_t data[28] = {0}, residual_data[4];
        BinaryMatrix h = {1+trial%4, 1+trial%7, data};
        if (trial == 1) { h.rows=3; h.cols=3; data[0]=data[4]=data[8]=1; } /* K=0. */
        else if (trial % 8) {
            for (size_t j = 0; j < h.rows*h.cols; ++j) { seed=1664525u*seed+1013904223u; data[j]=(uint8_t)(seed>>31); }
        } /* Every eighth trial has m=0. */
        BinarySystematic form = {0};
        CHECK(binary_matrix_systematize(&h, &form) == GF2_OK);
        unsigned target = syndrome(&h, 37);
        for (size_t i = 0; i < h.rows; ++i) residual_data[i]=(uint8_t)((target>>i)&1u);
        BitVector residual = {h.rows, residual_data};
        double weights[7];
        for (size_t j = 0; j < h.cols; ++j) weights[j] = trial%3 ? (double)((3*j+trial)%5) : 0;
        for (size_t delta = 0; delta <= form.rank; ++delta) verify(&h, &form, weights, &residual, delta);
        binary_systematic_free(&form);
    }
    printf("Local constraint checks passed (%u checks, %u partitions; Eq.9/10 exhaustive equivalence, all delta boundaries)\n", checks, partitions);
    return EXIT_SUCCESS;
}
