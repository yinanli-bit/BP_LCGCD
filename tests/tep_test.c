#include "bplcgcd/tep.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static unsigned checks, cases;
#define CHECK(condition) do { \
    if (!(condition)) { \
        fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, #condition); \
        exit(EXIT_FAILURE); \
    } \
    ++checks; \
} while (0)

typedef struct { unsigned mask, sorted_mask; double cost; } OraclePattern;

static int compare(const void *left, const void *right)
{
    const OraclePattern *a = left, *b = right;
    if (a->cost != b->cost) return a->cost < b->cost ? -1 : 1;
    return (a->sorted_mask > b->sorted_mask) - (a->sorted_mask < b->sorted_mask);
}

static unsigned mask_of(const BitVector *v)
{
    unsigned mask = 0;
    for (size_t j = 0; j < v->length; ++j) {
        CHECK(v->data[j] <= 1);
        mask |= (unsigned)v->data[j] << j;
    }
    return mask;
}

static void verify(const double *weights, size_t length, unsigned prefix)
{
    OraclePattern oracle[1024];
    unsigned rank[10], seen[1024] = {0}, total = 1u << length;
    double original[10];
    CHECK(length <= 10);
    if (length) memcpy(original, weights, length * sizeof(double));
    /* Brute-force every subset; determine sorted positions by counting lesser
     * (weight,index) pairs, independently of the production sort and FPT. */
    for (size_t j = 0; j < length; ++j) {
        rank[j] = 0;
        for (size_t i = 0; i < length; ++i) {
            if (weights[i] < weights[j] || (weights[i] == weights[j] && i < j)) ++rank[j];
        }
    }
    for (unsigned mask = 0; mask < total; ++mask) {
        oracle[mask] = (OraclePattern){.mask = mask};
        for (size_t r = 0; r < length; ++r) {
            for (size_t j = 0; j < length; ++j) {
                if (rank[j] == r && (mask & (1u << j))) {
                    oracle[mask].cost += weights[j];
                    oracle[mask].sorted_mask |= 1u << r;
                }
            }
        }
    }
    qsort(oracle, total, sizeof(*oracle), compare);
    TepGenerator g = {0};
    BitVector pattern = {0};
    CHECK(tep_generator_init(&g, weights, length) == TEP_OK);
    if (length) CHECK(bit_vector_init(&pattern, length) == GF2_OK);
    unsigned count = prefix && prefix < total ? prefix : total;
    double cost = -1, previous = -1;
    for (unsigned i = 0; i < count; ++i) {
        CHECK(tep_generator_next(&g, &pattern, &cost) == TEP_OK);
        unsigned mask = mask_of(&pattern);
        CHECK(mask == oracle[i].mask && cost == oracle[i].cost);
        CHECK(cost >= previous && !seen[mask]);
        if (!i) CHECK(mask == 0 && cost == 0);
        seen[mask] = 1;
        previous = cost;
    }
    if (count == total) {
        for (unsigned i = 0; i < total; ++i) CHECK(seen[i]);
        unsigned last = mask_of(&pattern);
        CHECK(tep_generator_next(&g, &pattern, &cost) == TEP_DONE);
        CHECK(mask_of(&pattern) == last && cost == previous);
        CHECK(tep_generator_next(&g, &pattern, &cost) == TEP_DONE);
    }
    if (length) CHECK(memcmp(weights, original, length * sizeof(double)) == 0);
    bit_vector_free(&pattern);
    tep_generator_free(&g); /* Also tests releasing an unconsumed frontier. */
    tep_generator_free(&g);
    CHECK(!g.nodes && !g.order && !g.size && !g.length);
    ++cases;
}

static void test_errors(void)
{
    TepGenerator g = {0};
    double weight = 1, cost = 9;
    CHECK(tep_generator_init(NULL, &weight, 1) == TEP_INVALID_ARGUMENT);
    CHECK(tep_generator_init(&g, NULL, 1) == TEP_INVALID_ARGUMENT);
    CHECK(tep_generator_init(&g, &weight, SIZE_MAX) == TEP_INVALID_ARGUMENT);
    const double bad[] = {-1, NAN, INFINITY};
    for (size_t i = 0; i < 3; ++i) {
        CHECK(tep_generator_init(&g, &bad[i], 1) == TEP_INVALID_ARGUMENT);
    }
    const double overflow[] = {DBL_MAX, DBL_MAX};
    CHECK(tep_generator_init(&g, overflow, 2) == TEP_NUMERIC_ERROR);
    CHECK(!g.nodes && !g.order);
    CHECK(tep_generator_init(&g, &weight, 1) == TEP_OK);
    CHECK(tep_generator_peek(NULL, &cost) == TEP_INVALID_ARGUMENT);
    CHECK(tep_generator_peek(&g, NULL) == TEP_INVALID_ARGUMENT);
    CHECK(tep_generator_peek(&g, &cost) == TEP_OK && cost == 0);
    cost = 9;
    BitVector pattern = {0};
    CHECK(bit_vector_init(&pattern, 1) == GF2_OK);
    pattern.data[0] = 1;
    CHECK(tep_generator_next(NULL, &pattern, &cost) == TEP_INVALID_ARGUMENT);
    CHECK(tep_generator_next(&g, NULL, &cost) == TEP_INVALID_ARGUMENT);
    CHECK(tep_generator_next(&g, &pattern, NULL) == TEP_INVALID_ARGUMENT);
    pattern.length = 2;
    CHECK(tep_generator_next(&g, &pattern, &cost) == TEP_DIMENSION_MISMATCH);
    CHECK(pattern.data[0] == 1 && cost == 9);
    pattern.length = 1;
    CHECK(tep_generator_next(&g, &pattern, &cost) == TEP_OK);
    CHECK(pattern.data[0] == 0 && cost == 0); /* Errors did not advance. */
    CHECK(g.nodes_expanded == 1);
    CHECK(tep_generator_peek(&g, &cost) == TEP_OK && cost == 1);
    CHECK(tep_generator_next(&g, &pattern, &cost) == TEP_OK && pattern.data[0] == 1);
    CHECK(g.nodes_expanded == 1);
    CHECK(tep_generator_peek(&g, &cost) == TEP_DONE && cost == 1);
    bit_vector_free(&pattern);
    tep_generator_free(&g);
    tep_generator_free(NULL);
}

static void test_wide_prefix(void)
{
    double weights[105], cost = -1, previous = -1;
    for (size_t j = 0; j < 105; ++j) weights[j] = 100 + (double)j;
    weights[104] = 0.25; weights[0] = 0.5; weights[64] = 0.75;
    TepGenerator g = {0};
    BitVector pattern = {0};
    CHECK(tep_generator_init(&g, weights, 105) == TEP_OK);
    CHECK(bit_vector_init(&pattern, 105) == GF2_OK);
    for (unsigned i = 0; i < 20; ++i) {
        CHECK(tep_generator_next(&g, &pattern, &cost) == TEP_OK);
        CHECK(cost >= previous);
        double actual = 0;
        for (size_t j = 0; j < 105; ++j) actual += pattern.data[j] * weights[j];
        CHECK(actual == cost);
        if (i == 0) CHECK(cost == 0);
        if (i == 1) CHECK(cost == 0.25 && pattern.data[104] == 1);
        if (i == 2) CHECK(cost == 0.5 && pattern.data[0] == 1);
        if (i == 3) CHECK(cost == 0.75 && pattern.data[0] == 1 && pattern.data[104] == 1);
        if (i == 4) CHECK(cost == 0.75 && pattern.data[64] == 1);
        previous = cost;
    }
    bit_vector_free(&pattern);
    tep_generator_free(&g);
    ++cases;
}

int main(void)
{
    test_errors();
    verify(NULL, 0, 0);
    double weights[10];
    uint32_t seed = 12345;
    for (size_t k = 1; k <= 10; ++k) {
        for (size_t j = 0; j < k; ++j) weights[j] = (double)(k - j);
        verify(weights, k, 0); /* Unsorted distinct weights. */
        for (size_t j = 0; j < k; ++j) weights[j] = 1;
        verify(weights, k, 0);
        for (size_t j = 0; j < k; ++j) weights[j] = 0;
        verify(weights, k, 0);
        for (size_t j = 0; j < k; ++j) {
            seed = 1664525u * seed + 1013904223u;
            weights[j] = (double)((seed >> 24) % 13) / 4;
        }
        verify(weights, k, 0); /* Binary fractions: exact cost comparisons. */
        verify(weights, k, 17);
    }
    const double example[] = {3, 1, 2};
    verify(example, 3, 0);
    const double scales[] = {1e16, 1, 1, 0.25, 0};
    verify(scales, 5, 0); /* Avoid cancellation in incremental cost updates. */
    const double fractional[] = {0.3, 0.1, 0.2, 0.1};
    verify(fractional, 4, 0);
    test_wide_prefix();
    printf("TEP checks passed (%u checks, %u cases; exhaustive K<=10, K=105 prefix)\n", checks, cases);
    return EXIT_SUCCESS;
}
