#include "bplcgcd/decoder.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static unsigned checks;
#define CHECK(condition) do { \
    if (!(condition)) { \
        fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, #condition); \
        exit(EXIT_FAILURE); \
    } \
    ++checks; \
} while (0)

/* Independent integer oracle for the three rows (column zero is unused). */
static unsigned syndrome_bits(unsigned x)
{
    unsigned b1 = (x >> 1) & 1u, b2 = (x >> 2) & 1u, b3 = (x >> 3) & 1u;
    return ((b1 + b2) % 2) | (((b2 + b3) % 2) << 1) | (((b1 + b3) % 2) << 2);
}

static void set_bits(BitVector *v, unsigned bits)
{
    for (size_t j = 0; j < v->length; ++j) v->data[j] = (uint8_t)((bits >> j) & 1u);
}

int main(void)
{
    BinaryMatrix h = {0};
    BitVector reference = {0}, residual = {0}, target = {0}, hint = {0}, pattern = {0}, scratch = {0};
    DecodeResult result = {0};
    CHECK(binary_matrix_init(&h, 3, 4) == GF2_OK);
    const uint8_t entries[] = {0,1,1,0, 0,0,1,1, 0,1,0,1};
    memcpy(h.data, entries, sizeof(entries));
    CHECK(bit_vector_init(&reference, 4) == GF2_OK);
    CHECK(bit_vector_init(&residual, 3) == GF2_OK);
    CHECK(bit_vector_init(&target, 3) == GF2_OK);
    CHECK(bit_vector_init(&hint, 4) == GF2_OK);
    CHECK(bit_vector_init(&pattern, 4) == GF2_OK);
    CHECK(bit_vector_init(&scratch, 3) == GF2_OK);
    CHECK(bit_vector_init(&result.estimate, 4) == GF2_OK);
    double llr[] = {-0.0, -2, 3, 0.0}, weights[] = {5, 2, 3, 7};
    double flip_weights[4];
    DecodeInput input = {.mode=DECODE_CLASSICAL, .h=&h, .length=4, .llr=llr};
    CHECK(decode_prepare(&input, &reference, &residual, flip_weights) == DECODE_SUCCESS);
    CHECK(reference.data[0] == 0 && reference.data[1] == 1 && reference.data[3] == 0);
    CHECK(residual.data[0] == 1 && residual.data[1] == 0 && residual.data[2] == 1);
    CHECK(flip_weights[0] == 0 && flip_weights[1] == 2 && flip_weights[2] == 3 && flip_weights[3] == 0);
    for (unsigned x = 0; x < 16; ++x) {
        set_bits(&pattern, x);
        DecodeStatus status = decode_check_candidate(&input, flip_weights, &pattern, &result, &scratch);
        CHECK(status == (syndrome_bits(x ^ 2u) == 0 ? DECODE_SUCCESS : DECODE_INVALID_CANDIDATE));
        CHECK(!result.optimality_proven);
        CHECK(result.candidate_valid == (status == DECODE_SUCCESS));
        if (result.candidate_valid) {
            CHECK(result.cost == 2 * ((x >> 1) & 1u) + 3 * ((x >> 2) & 1u));
            for (size_t j = 0; j < 4; ++j) CHECK(result.estimate.data[j] == (((x ^ 2u) >> j) & 1u));
        }
    }
    set_bits(&pattern, 2);
    CHECK(decode_check_candidate(&input, flip_weights, &pattern, &result, &scratch) == DECODE_SUCCESS);
    CHECK(result.cost == 2);
    result.stats.candidates_popped = 9;
    result.stats.candidates_evaluated = 4;
    result.stats.nodes_expanded = 6;
    result.stats.batches = 2;
    result.optimality_proven = true; /* Simulate a previously certified result. */
    decode_result_mark_truncated(&result);
    CHECK(result.status == DECODE_TRUNCATED && result.candidate_valid && !result.optimality_proven);
    CHECK(result.cost == 2 && result.stats.candidates_popped == 9 && result.stats.batches == 2);
    uint8_t *buffer = result.estimate.data;
    decode_result_reset(&result);
    CHECK(result.estimate.data == buffer && result.estimate.length == 4);
    CHECK(result.status == DECODE_NOT_RUN && !result.candidate_valid && !result.optimality_proven);
    CHECK(!result.stats.candidates_evaluated && !result.stats.nodes_expanded);
    decode_result_mark_truncated(&result);
    CHECK(result.status == DECODE_TRUNCATED && !result.candidate_valid);

    input = (DecodeInput){.mode=DECODE_SYNDROME, .h=&h, .length=4,
                          .weights=weights, .s_target=&target, .z_hint=&hint};
    for (unsigned t = 0; t < 8; ++t) {
        set_bits(&target, t);
        for (unsigned z = 0; z < 16; ++z) {
            set_bits(&hint, z);
            CHECK(decode_prepare(&input, &reference, &residual, flip_weights) == DECODE_SUCCESS);
            CHECK(memcmp(reference.data, hint.data, 4) == 0);
            for (size_t i = 0; i < 3; ++i) CHECK(residual.data[i] == (((t ^ syndrome_bits(z)) >> i) & 1u));
            for (unsigned f = 0; f < 16; ++f) {
                set_bits(&pattern, f);
                DecodeStatus status = decode_check_candidate(&input, flip_weights, &pattern, &result, &scratch);
                CHECK(status == (syndrome_bits(z ^ f) == t ? DECODE_SUCCESS : DECODE_INVALID_CANDIDATE));
                CHECK(result.candidate_valid == (status == DECODE_SUCCESS) && !result.optimality_proven);
                if (result.candidate_valid) {
                    double cost = 0;
                    for (size_t j = 0; j < 4; ++j) {
                        unsigned bit = ((z ^ f) >> j) & 1u;
                        CHECK(result.estimate.data[j] == bit);
                        cost += ((f >> j) & 1u) * flip_weights[j];
                    }
                    CHECK(result.cost == cost);
                }
            }
        }
    }
    set_bits(&target, 5); set_bits(&hint, 5); set_bits(&pattern, 7);
    CHECK(decode_prepare(&input, &reference, &residual, flip_weights) == DECODE_SUCCESS);
    CHECK(residual.data[0] == 0 && residual.data[1] == 1 && residual.data[2] == 1);
    CHECK(decode_check_candidate(&input, flip_weights, &pattern, &result, &scratch) == DECODE_SUCCESS);
    CHECK(result.cost == 10 && result.estimate.data[1] == 1); /* Score f, not final e_hat (cost 2). */
    input.z_hint = NULL;
    CHECK(decode_prepare(&input, &reference, &residual, flip_weights) == DECODE_SUCCESS);
    CHECK(reference.data[0] == 0 && reference.data[1] == 0 && reference.data[2] == 0 && reference.data[3] == 0);
    CHECK(memcmp(residual.data, target.data, 3) == 0);
    CHECK(memcmp(flip_weights, weights, sizeof(weights)) == 0);
    set_bits(&pattern, 2);
    CHECK(decode_check_candidate(&input, flip_weights, &pattern, &result, &scratch) == DECODE_SUCCESS);
    CHECK(result.cost == 2 && memcmp(result.estimate.data, pattern.data, 4) == 0);
    weights[1] = 0; set_bits(&pattern, 2);
    CHECK(decode_prepare(&input, &reference, &residual, flip_weights) == DECODE_SUCCESS);
    CHECK(decode_check_candidate(&input, flip_weights, &pattern, &result, &scratch) == DECODE_SUCCESS && result.cost == 0);
    weights[1] = 2;
    weights[3] = -1;
    CHECK(decode_prepare(&input, &reference, &residual, flip_weights) == DECODE_INVALID_ARGUMENT);
    weights[3] = NAN;
    CHECK(decode_prepare(&input, &reference, &residual, flip_weights) == DECODE_INVALID_ARGUMENT);
    weights[3] = 7;
    input.weights = NULL;
    CHECK(decode_prepare(&input, &reference, &residual, flip_weights) == DECODE_INVALID_ARGUMENT);
    input.weights = weights; input.s_target = NULL;
    CHECK(decode_prepare(&input, &reference, &residual, flip_weights) == DECODE_INVALID_ARGUMENT);
    input.s_target = &target; target.length = 2;
    CHECK(decode_prepare(&input, &reference, &residual, flip_weights) == DECODE_DIMENSION_MISMATCH);
    target.length = 3; hint.length = 3; input.z_hint = &hint;
    CHECK(decode_prepare(&input, &reference, &residual, flip_weights) == DECODE_DIMENSION_MISMATCH);
    hint.length = 4; input.z_hint = NULL;
    input.length = 3;
    CHECK(decode_prepare(&input, &reference, &residual, flip_weights) == DECODE_DIMENSION_MISMATCH);
    input.length = 4; input.llr = llr;
    CHECK(decode_prepare(&input, &reference, &residual, flip_weights) == DECODE_INVALID_ARGUMENT);
    input.llr = NULL;
    CHECK(decode_prepare(&input, &reference, &residual, NULL) == DECODE_INVALID_ARGUMENT);

    /* BP LLR is a source adapter, not an implemented BP decoder. */
    double bp_llr[] = {2, -5, 1, -7};
    input = (DecodeInput){.mode=DECODE_BP_SYNDROME, .h=&h, .length=4,
                          .llr=bp_llr, .s_target=&target};
    CHECK(decode_prepare(&input, &reference, &residual, flip_weights) == DECODE_SUCCESS);
    set_bits(&hint, 10);
    CHECK(memcmp(reference.data, hint.data, 4) == 0);
    CHECK(residual.data[0] == 0 && residual.data[1] == 1 && residual.data[2] == 1);
    double constant = 0;
    for (size_t j = 0; j < 4; ++j) {
        CHECK(flip_weights[j] == fabs(bp_llr[j]));
        if (hint.data[j]) constant += bp_llr[j];
    }
    for (unsigned f = 0; f < 16; ++f) {
        double signed_score = 0, flip_cost = 0;
        set_bits(&pattern, f);
        for (size_t j = 0; j < 4; ++j) {
            signed_score += (hint.data[j] ^ pattern.data[j]) * bp_llr[j];
            flip_cost += pattern.data[j] * flip_weights[j];
        }
        CHECK(signed_score == constant + flip_cost);
        DecodeStatus status = decode_check_candidate(&input, flip_weights, &pattern, &result, &scratch);
        CHECK(status == (syndrome_bits(10u ^ f) == 5 ? DECODE_SUCCESS : DECODE_INVALID_CANDIDATE));
        if (status == DECODE_SUCCESS) CHECK(result.cost == flip_cost);
    }
    set_bits(&pattern, 6); /* e_hat=(0,0,1,1): correct W=6, wrong absolute score=8. */
    CHECK(decode_check_candidate(&input, flip_weights, &pattern, &result, &scratch) == DECODE_SUCCESS);
    CHECK(result.cost == 6 && result.estimate.data[2] == 1 && result.estimate.data[3] == 1);
    double first_cost = result.cost, first_absolute_score = 0;
    for (size_t j = 0; j < 4; ++j) first_absolute_score += result.estimate.data[j] * flip_weights[j];
    set_bits(&pattern, 8); /* e_hat=(0,1,0,0): correct W=7, wrong absolute score=5. */
    CHECK(decode_check_candidate(&input, flip_weights, &pattern, &result, &scratch) == DECODE_SUCCESS);
    CHECK(result.cost == 7 && result.estimate.data[1] == 1); /* Wrong scoring reverses the ranking. */
    double second_absolute_score = 0;
    for (size_t j = 0; j < 4; ++j) second_absolute_score += result.estimate.data[j] * flip_weights[j];
    CHECK(first_absolute_score == 8 && second_absolute_score == 5);
    CHECK(first_cost < result.cost && first_absolute_score > second_absolute_score);
    input = (DecodeInput){.mode=DECODE_SYNDROME, .h=&h, .length=4,
                          .weights=flip_weights, .s_target=&target, .z_hint=&hint};
    double explicit_weights[4];
    CHECK(decode_prepare(&input, &reference, &residual, explicit_weights) == DECODE_SUCCESS);
    CHECK(memcmp(explicit_weights, flip_weights, sizeof(flip_weights)) == 0);
    CHECK(decode_check_candidate(&input, explicit_weights, &pattern, &result, &scratch) == DECODE_SUCCESS && result.cost == 7);

    input = (DecodeInput){.mode=DECODE_CLASSICAL, .h=&h, .length=4, .llr=llr};
    llr[3] = INFINITY;
    CHECK(decode_prepare(&input, &reference, &residual, flip_weights) == DECODE_INVALID_ARGUMENT);
    llr[3] = 0; input.weights = weights;
    CHECK(decode_prepare(&input, &reference, &residual, flip_weights) == DECODE_INVALID_ARGUMENT);
    input.weights = NULL; input.mode = (DecodeMode)99;
    CHECK(decode_prepare(&input, &reference, &residual, flip_weights) == DECODE_INVALID_ARGUMENT);
    input.mode = DECODE_CLASSICAL;
    CHECK(decode_prepare(NULL, &reference, &residual, flip_weights) == DECODE_INVALID_ARGUMENT);
    input.llr = NULL;
    CHECK(decode_prepare(&input, &reference, &residual, flip_weights) == DECODE_INVALID_ARGUMENT);
    input.llr = llr;
    CHECK(decode_prepare(&input, &reference, &residual, flip_weights) == DECODE_SUCCESS);
    set_bits(&pattern, 2);
    CHECK(decode_check_candidate(&input, flip_weights, &pattern, &result, &scratch) == DECODE_SUCCESS);
    set_bits(&pattern, 0);
    CHECK(decode_check_candidate(&input, flip_weights, &pattern, &result, &scratch) == DECODE_INVALID_CANDIDATE);
    CHECK(!result.candidate_valid && !result.optimality_proven);
    CHECK(decode_check_candidate(&input, flip_weights, &pattern, NULL, &scratch) == DECODE_INVALID_ARGUMENT);
    pattern.data[3] = 0; result.estimate.length = 3;
    CHECK(decode_check_candidate(&input, flip_weights, &pattern, &result, &scratch) == DECODE_DIMENSION_MISMATCH);
    result.estimate.length = 4;
    llr[0] = DBL_MAX; llr[1] = -DBL_MAX; set_bits(&pattern, 3);
    CHECK(decode_prepare(&input, &reference, &residual, flip_weights) == DECODE_SUCCESS);
    CHECK(decode_check_candidate(&input, flip_weights, &pattern, &result, &scratch) == DECODE_NUMERIC_ERROR);
    CHECK(!result.candidate_valid);
    llr[0] = -0.0; llr[1] = -2;
    llr[1] = llr[2] = 0; set_bits(&pattern, 0);
    CHECK(decode_prepare(&input, &reference, &residual, flip_weights) == DECODE_SUCCESS);
    CHECK(decode_check_candidate(&input, flip_weights, &pattern, &result, &scratch) == DECODE_SUCCESS);
    CHECK(result.cost == 0 && result.candidate_valid && !result.optimality_proven);
    llr[1] = -2; llr[2] = 3;
    CHECK(decode_prepare(&input, &reference, &residual, flip_weights) == DECODE_SUCCESS);
    const char *names[] = {"gcd", "lcgcd", "bp", "nms", "bp-lcgcd", "bp-osd"};
    for (size_t i = 0; i < 6; ++i) {
        DecodeConfig config = {.decoder_name=names[i], .candidate_limit=1};
        CHECK(decoder_is_known(names[i]));
        set_bits(&pattern, 2);
        CHECK(decode_check_candidate(&input, flip_weights, &pattern, &result, &scratch) == DECODE_SUCCESS);
        result.stats.candidates_evaluated = 4;
        DecodeStatus expected = i == 0 ? DECODE_INVALID_ARGUMENT : DECODE_NOT_IMPLEMENTED; /* GCD needs S03 form. */
        CHECK(decoder_decode(&config, &input, &result) == expected);
        CHECK(result.status == expected && !result.candidate_valid && !result.optimality_proven);
        CHECK(result.stats.candidates_evaluated == 0 && result.estimate.data == buffer);
    }
    CHECK(!decoder_is_known(NULL) && !decoder_is_known("GCD") && !decoder_is_known("unknown"));
    DecodeConfig unknown = {.decoder_name="unknown"};
    CHECK(decoder_decode(&unknown, &input, &result) == DECODE_INVALID_ARGUMENT);
    CHECK(decoder_decode(NULL, &input, &result) == DECODE_INVALID_ARGUMENT);
    CHECK(decoder_decode(&unknown, &input, NULL) == DECODE_INVALID_ARGUMENT);
    decode_result_reset(NULL); decode_result_mark_truncated(NULL);
    bit_vector_free(&result.estimate);
    bit_vector_free(&reference); bit_vector_free(&residual); bit_vector_free(&target);
    bit_vector_free(&hint); bit_vector_free(&pattern); bit_vector_free(&scratch); binary_matrix_free(&h);
    printf("Decoder interface checks passed (%u checks; unified flip cost in all three modes)\n", checks);
    return EXIT_SUCCESS;
}
