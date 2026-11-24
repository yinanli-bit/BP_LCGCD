#include "bplcgcd/decoder.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

static unsigned checks, frames, proven, truncated, inconsistent;
static uint32_t seed = 78901;
#define CHECK(condition) do { \
    if (!(condition)) { \
        fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, #condition); \
        exit(EXIT_FAILURE); \
    } \
    ++checks; \
} while (0)

static unsigned random_value(void)
{
    seed = 1664525u * seed + 1013904223u;
    return seed >> 16;
}

static unsigned syndrome(const BinaryMatrix *h, unsigned bits)
{
    unsigned value = 0;
    for (size_t i = 0; i < h->rows; ++i) {
        unsigned sum = 0;
        for (size_t j = 0; j < h->cols; ++j) sum += h->data[i * h->cols + j] * ((bits >> j) & 1u);
        value |= (sum % 2) << i;
    }
    return value;
}

static unsigned bits_of(const BitVector *v)
{
    unsigned bits = 0;
    for (size_t j = 0; j < v->length; ++j) {
        CHECK(v->data[j] <= 1);
        bits |= (unsigned)v->data[j] << j;
    }
    return bits;
}

static void verify(const DecodeInput *input, uint64_t limit)
{
    const BinaryMatrix *h = input->h;
    unsigned reference = 0, target = input->mode == DECODE_CLASSICAL ? 0 : bits_of(input->s_target);
    double weights[8];
    for (size_t j = 0; j < h->cols; ++j) {
        unsigned bit = input->mode == DECODE_SYNDROME ?
                       (input->z_hint ? input->z_hint->data[j] : 0) : (unsigned)(input->llr[j] < 0);
        reference |= bit << j;
        weights[j] = input->mode == DECODE_SYNDROME ? input->weights[j] : fabs(input->llr[j]);
    }
    /* Independent ML oracle: enumerate FULL original-coordinate flips, without
     * S03 transforms, TEP ordering, basis completion or production scoring. */
    double optimum = INFINITY;
    for (unsigned f = 0; f < (1u << h->cols); ++f) {
        if (syndrome(h, reference ^ f) != target) continue;
        double cost = 0;
        for (size_t j = 0; j < h->cols; ++j) cost += ((f >> j) & 1u) * weights[j];
        if (cost < optimum) optimum = cost;
    }
    DecodeConfig config = {.decoder_name="gcd", .candidate_limit=limit};
    DecodeResult result = {0};
    CHECK(bit_vector_init(&result.estimate, h->cols) == GF2_OK);
    DecodeStatus status = decoder_decode(&config, input, &result);
    if (!isfinite(optimum)) {
        ++inconsistent;
        CHECK(status == DECODE_INCONSISTENT_SYNDROME);
        CHECK(!result.candidate_valid && !result.optimality_proven);
        CHECK(!result.stats.candidates_popped && !result.stats.candidates_evaluated);
        CHECK(result.stop_reason == DECODE_STOP_NONE);
    } else {
        CHECK(status == DECODE_SUCCESS || status == DECODE_TRUNCATED);
        CHECK(result.status == status && result.candidate_valid);
        unsigned estimate = bits_of(&result.estimate), f = reference ^ estimate;
        CHECK(syndrome(h, estimate) == target);
        double cost = 0;
        for (size_t j = 0; j < h->cols; ++j) cost += ((f >> j) & 1u) * weights[j];
        CHECK(result.cost == cost && cost >= optimum);
        CHECK(result.stats.candidates_popped == result.stats.candidates_evaluated);
        CHECK(result.stats.nodes_expanded <= result.stats.candidates_popped && result.stats.batches == 0);
        CHECK(!limit || result.stats.candidates_evaluated <= limit);
        if (!limit) CHECK(status == DECODE_SUCCESS);
        if (status == DECODE_SUCCESS) {
            ++proven;
            CHECK(result.optimality_proven && cost == optimum);
            CHECK(result.stop_reason == DECODE_STOP_BOUND || result.stop_reason == DECODE_STOP_EXHAUSTED);
        } else {
            ++truncated;
            CHECK(!result.optimality_proven && result.stop_reason == DECODE_STOP_BUDGET);
        }
    }
    bit_vector_free(&result.estimate);
    ++frames;
}

static void exercise(BinaryMatrix *h, unsigned trial)
{
    BinarySystematic form = {0};
    BitVector hint = {0}, target = {0};
    CHECK(binary_matrix_systematize(h, &form) == GF2_OK); /* Offline, once per H. */
    CHECK(bit_vector_init(&hint, h->cols) == GF2_OK);
    CHECK(bit_vector_init(&target, h->rows) == GF2_OK);
    double llr[8], weights[8];
    for (size_t j = 0; j < h->cols; ++j) {
        weights[j] = (double)(random_value() % 17) / 4;
        llr[j] = random_value() & 1u ? -weights[j] : weights[j];
        hint.data[j] = (uint8_t)(random_value() & 1u);
    }
    unsigned target_bits = trial & 1u ? random_value() & ((1u << h->rows) - 1) :
                           syndrome(h, random_value() & ((1u << h->cols) - 1));
    for (size_t i = 0; i < h->rows; ++i) target.data[i] = (uint8_t)((target_bits >> i) & 1u);
    for (unsigned mode = 0; mode < 3; ++mode) {
        DecodeInput input = {.mode=(DecodeMode)mode, .h=h, .length=h->cols, .systematic=&form};
        if (mode == DECODE_SYNDROME) {
            input.weights = weights;
            if (trial & 2u) input.z_hint = &hint;
        } else input.llr = llr;
        if (mode != DECODE_CLASSICAL) input.s_target = &target;
        verify(&input, 0);
        verify(&input, 1);
    }
    binary_systematic_free(&form);
    bit_vector_free(&hint);
    bit_vector_free(&target);
}

static void test_budget_and_bound(void)
{
    BinaryMatrix h = {0};
    BinarySystematic form = {0};
    BitVector target = {0};
    DecodeResult result = {0};
    CHECK(binary_matrix_init(&h, 2, 4) == GF2_OK);
    const uint8_t entries[] = {1,0,1,0, 0,1,0,1};
    memcpy(h.data, entries, sizeof(entries));
    CHECK(binary_matrix_systematize(&h, &form) == GF2_OK);
    CHECK(bit_vector_init(&target, 2) == GF2_OK);
    CHECK(bit_vector_init(&result.estimate, 4) == GF2_OK);
    target.data[0] = target.data[1] = 1;
    double weights[] = {10,10,1,1};
    DecodeInput input = {.mode=DECODE_SYNDROME, .h=&h, .length=4,
                          .weights=weights, .s_target=&target, .systematic=&form};
    DecodeConfig config = {.decoder_name="gcd", .candidate_limit=1};
    CHECK(decoder_decode(&config, &input, &result) == DECODE_TRUNCATED);
    CHECK(result.candidate_valid && !result.optimality_proven && result.cost == 20);
    CHECK(result.stop_reason == DECODE_STOP_BUDGET);
    CHECK(result.stats.candidates_evaluated == 1 && result.stats.candidates_popped == 1);
    config.candidate_limit = 2;
    CHECK(decoder_decode(&config, &input, &result) == DECODE_TRUNCATED && result.cost == 11);
    config.candidate_limit = 0;
    CHECK(decoder_decode(&config, &input, &result) == DECODE_SUCCESS);
    CHECK(result.candidate_valid && result.optimality_proven && result.cost == 2);
    CHECK(bits_of(&result.estimate) == 12);
    CHECK(result.stats.candidates_evaluated == 4 && result.stats.nodes_expanded == 2);
    CHECK(result.stop_reason == DECODE_STOP_EXHAUSTED);
    config.candidate_limit = 4; /* Exhaustion proves optimal even at the budget. */
    CHECK(decoder_decode(&config, &input, &result) == DECODE_SUCCESS && result.optimality_proven);
    target.data[1] = 0; weights[0] = 0.5; config.candidate_limit = 1;
    CHECK(decoder_decode(&config, &input, &result) == DECODE_SUCCESS);
    CHECK(result.optimality_proven && result.cost == 0.5 && result.stats.candidates_evaluated == 1);
    CHECK(result.stop_reason == DECODE_STOP_BOUND);
    for (size_t j = 0; j < 4; ++j) weights[j] = 0;
    CHECK(decoder_decode(&config, &input, &result) == DECODE_SUCCESS && result.cost == 0);
    unsigned first = bits_of(&result.estimate);
    CHECK(decoder_decode(&config, &input, &result) == DECODE_SUCCESS && bits_of(&result.estimate) == first);
    bit_vector_free(&result.estimate); bit_vector_free(&target);
    binary_systematic_free(&form); binary_matrix_free(&h);
}

int main(void)
{
    test_budget_and_bound();
    for (unsigned trial = 0; trial < 64; ++trial) {
        BinaryMatrix h = {0};
        CHECK(binary_matrix_init(&h, 1 + trial % 5, 1 + trial % 8) == GF2_OK);
        for (size_t i = 0; i < h.rows * h.cols; ++i) h.data[i] = (uint8_t)(random_value() & 1u);
        if (trial % 8 == 0) memset(h.data, 0, h.rows * h.cols); /* rank=0. */
        exercise(&h, trial);
        binary_matrix_free(&h);
    }
    BinaryMatrix h = {0};
    CHECK(binary_matrix_init(&h, 4, 7) == GF2_OK);
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 7; ++j) h.data[i * 7 + j] = (uint8_t)(((j + 1) >> i) & 1u);
    }
    for (size_t j = 0; j < 7; ++j) h.data[21 + j] = h.data[j] ^ h.data[7 + j];
    exercise(&h, 0); /* Hamming [7,4] plus a redundant row, with column swaps. */
    exercise(&h, 1); /* Also compare possibly inconsistent syndromes. */
    binary_matrix_free(&h);
    CHECK(binary_matrix_init(&h, 2, 2) == GF2_OK);
    h.data[0] = h.data[3] = 1;
    exercise(&h, 0); /* Full column rank, K=0: only the empty partial pattern. */
    binary_matrix_free(&h);
    CHECK(proven && truncated && inconsistent);
    printf("GCD checks passed (%u checks, %u oracle frames; proven=%u, truncated=%u, inconsistent=%u)\n",
           checks, frames, proven, truncated, inconsistent);
    return EXIT_SUCCESS;
}
