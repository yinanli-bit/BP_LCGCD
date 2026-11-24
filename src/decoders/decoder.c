#include "bplcgcd/decoder.h"
#include "bplcgcd/tep.h"

#include <inttypes.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static DecodeStatus vector_status(const BitVector *v, size_t length)
{
    if (!v || !v->data) return DECODE_INVALID_ARGUMENT;
    if (v->length != length) return DECODE_DIMENSION_MISMATCH;
    return DECODE_SUCCESS;
}

static DecodeStatus validate_input(const DecodeInput *input)
{
    if (!input || !input->h) return DECODE_INVALID_ARGUMENT;
    const BinaryMatrix *h = input->h;
    if (!h->data || !h->rows || !h->cols) {
        return DECODE_INVALID_ARGUMENT;
    }
    if (input->length != h->cols) return DECODE_DIMENSION_MISMATCH;
    const double *values;
    if (input->mode == DECODE_CLASSICAL) {
        if (!input->llr || input->weights || input->s_target || input->z_hint) {
            return DECODE_INVALID_ARGUMENT;
        }
        values = input->llr;
    } else if (input->mode == DECODE_SYNDROME || input->mode == DECODE_BP_SYNDROME) {
        DecodeStatus status = vector_status(input->s_target, h->rows);
        if (status != DECODE_SUCCESS) return status;
        if (input->mode == DECODE_SYNDROME) {
            if (input->llr || !input->weights) return DECODE_INVALID_ARGUMENT;
            if (input->z_hint) {
                status = vector_status(input->z_hint, h->cols);
                if (status != DECODE_SUCCESS) return status;
            }
            values = input->weights;
        } else {
            if (!input->llr || input->weights || input->z_hint) return DECODE_INVALID_ARGUMENT;
            values = input->llr;
        }
    } else {
        return DECODE_INVALID_ARGUMENT;
    }
    for (size_t j = 0; j < input->length; ++j) {
        if (!isfinite(values[j]) ||
            (input->mode == DECODE_SYNDROME && values[j] < 0)) {
            return DECODE_INVALID_ARGUMENT;
        }
    }
    return DECODE_SUCCESS;
}

//硬判决或0
static uint8_t reference_bit(const DecodeInput *input, size_t j)
{
    if (input->mode != DECODE_SYNDROME) return (uint8_t)(input->llr[j] < 0);
    return input->z_hint ? input->z_hint->data[j] : 0;
} 

//得到归一化搜索输入s_res和翻转代价flip_weights（是llr的绝对值或weights）
DecodeStatus decode_prepare(const DecodeInput *input, BitVector *reference,BitVector *s_res, double *flip_weights)
{
    DecodeStatus status = validate_input(input);
    if (status != DECODE_SUCCESS) return status;
    if (!reference || !s_res || !reference->data || !s_res->data || !flip_weights) {
        return DECODE_INVALID_ARGUMENT;
    }
    if (reference->length != input->length || s_res->length != input->h->rows) {
        return DECODE_DIMENSION_MISMATCH;
    }
    for (size_t j = 0; j < input->length; ++j) {
        reference->data[j] = reference_bit(input, j);
        flip_weights[j] = input->mode == DECODE_SYNDROME ? input->weights[j] : fabs(input->llr[j]); //物理错误先验 or 信道llr or BP后验llr
    }
    if (binary_matrix_multiply_vector(input->h, reference, s_res) != GF2_OK) {
        return DECODE_INVALID_ARGUMENT;
    }
    if (input->mode != DECODE_CLASSICAL) {
        for (size_t i = 0; i < s_res->length; ++i) {
            s_res->data[i] ^= input->s_target->data[i];
        }
    }
    return DECODE_SUCCESS;
}

void decode_result_reset(DecodeResult *result)
{
    if (result) {
        BitVector buffer = result->estimate;
        *result = (DecodeResult){.estimate = buffer};
    }
}// reset result 仅保留 estimate buffer，其他成员置为0

void decode_result_mark_truncated(DecodeResult *result)
{
    if (result) {
        result->status = DECODE_TRUNCATED;
        result->optimality_proven = false;
        result->stop_reason = DECODE_STOP_BUDGET;
    }
}

//根据候选 f，恢复估计estimate，计算代价 cost，检查候选是否有效
DecodeStatus decode_check_candidate(const DecodeInput *input,
                                    const double *flip_weights, const BitVector *f,
                                    DecodeResult *result, BitVector *syndrome_scratch)
{
    if (!result) return DECODE_INVALID_ARGUMENT;
    result->candidate_valid = result->optimality_proven = false;
    result->cost = 0;
    result->stop_reason = DECODE_STOP_NONE;
    if (!input || !input->h || !flip_weights || !f || !f->data ||
        !result->estimate.data || !syndrome_scratch || !syndrome_scratch->data) {
        return result->status = DECODE_INVALID_ARGUMENT;
    }
    if (f->length != input->length || result->estimate.length != input->length ||
        syndrome_scratch->length != input->h->rows) {
        return result->status = DECODE_DIMENSION_MISMATCH;
    }
    for (size_t j = 0; j < input->length; ++j) {
        result->estimate.data[j] = reference_bit(input, j) ^ f->data[j];
    }
    if (binary_matrix_multiply_vector(input->h, &result->estimate, syndrome_scratch) != GF2_OK) {
        return result->status = DECODE_INVALID_ARGUMENT;
    }
    for (size_t i = 0; i < input->h->rows; ++i) {
        uint8_t target = input->mode == DECODE_CLASSICAL ? 0 : input->s_target->data[i];
        if (syndrome_scratch->data[i] != target) return result->status = DECODE_INVALID_CANDIDATE;
    }
    //可以一边算一边判断。

    double cost = 0;
    for (size_t j = 0; j < input->length; ++j) {
        if (f->data[j]) cost += flip_weights[j];
    }
    if (!isfinite(cost)) return result->status = DECODE_NUMERIC_ERROR;
    result->cost = cost;
    result->candidate_valid = true;
    return result->status = DECODE_SUCCESS;
}

static void trace_bits(FILE *stream, const char *label, const uint8_t *bits, size_t length)
{
    fputs(label, stream);
    if (!length) fputc('-', stream);
    for (size_t j = 0; j < length; ++j) fputc('0' + bits[j], stream);
}

static DecodeStatus gcd_decode(const DecodeConfig *config, const DecodeInput *input,DecodeResult *result)
{
    if (!input || !input->h || !input->h->data || !result->estimate.data) {
        return result->status = DECODE_INVALID_ARGUMENT;
    }
    const BinarySystematic *form = input->systematic;
    if (!form || !form->reduced.data || !form->row_transform.data ||
        !form->sys_to_orig || !form->orig_to_sys) return result->status = DECODE_INVALID_ARGUMENT;
    size_t n = input->h->cols, rows = input->h->rows, m = form->rank;
    if (!n || !rows || n > SIZE_MAX / sizeof(double)) return result->status = DECODE_INVALID_ARGUMENT;
    if (result->estimate.length != n || form->reduced.rows != rows || form->reduced.cols != n ||
        form->row_transform.rows != rows || form->row_transform.cols != rows || m > n || m > rows) {
        return result->status = DECODE_DIMENSION_MISMATCH;
    }
    size_t k = n - m;
    BitVector reference = {0}, residual = {0}, transformed = {0}, system = {0};
    BitVector partial = {0}, original = {0}, scratch = {0};
    BitVector *buffers[] = {&reference, &residual, &transformed, &system, &partial, &original, &scratch};
    const size_t lengths[] = {n, rows, rows, n, k, n, rows};
    TepGenerator generator = {0};
    TepStatus tep_status = TEP_OK;
    DecodeStatus status = DECODE_ALLOCATION_FAILED;
    double *weights = malloc(n * sizeof(double));
    double *partial_weights = k ? malloc(k * sizeof(double)) : NULL;
    if (!weights || (k && !partial_weights)) goto finish;
    for (size_t i = 0; i < sizeof(buffers) / sizeof(*buffers); ++i) {
        if (lengths[i] && bit_vector_init(buffers[i], lengths[i]) != GF2_OK) goto finish;
    }
    status = decode_prepare(input, &reference, &residual, weights);
    if (status != DECODE_SUCCESS) goto finish;
    Gf2Status gf2_status = binary_systematic_transform_syndrome(form, &residual, &transformed);
    if (gf2_status != GF2_OK) {
        status = gf2_status == GF2_INCONSISTENT_SYNDROME ?
                 DECODE_INCONSISTENT_SYNDROME : DECODE_INVALID_ARGUMENT;
        goto finish;
    }
    for (size_t j = 0; j < k; ++j) partial_weights[j] = weights[form->sys_to_orig[m + j]];
    tep_status = tep_generator_init(&generator, partial_weights, k);
    

    if (tep_status != TEP_OK) goto tep_error;
    FILE *trace = config->trace;
    if (trace) {
        fprintf(trace, "trace N=%zu rank=%zu K=%zu sys_to_orig=", n, m, k);
        for (size_t j = 0; j < n; ++j) fprintf(trace, "%s%zu", j ? "," : "", form->sys_to_orig[j]);
        trace_bits(trace, "\ntrace reference(original)=", reference.data, n);
        trace_bits(trace, "\ntrace s_res(original_rows)=", residual.data, rows);
        trace_bits(trace, "\ntrace U_s_res=", transformed.data, rows);
        fputs("\ntrace flip_weights(original)=", trace);
        for (size_t j = 0; j < n; ++j) fprintf(trace, "%s%.17g", j ? "," : "", weights[j]);
        fputc('\n', trace);
        for (size_t i = 0; i < m; ++i) {
            fprintf(trace, "trace H_sys[%zu]=", i);
            trace_bits(trace, "", form->reduced.data + i * n, n);
            fputc('\n', trace);
        }
    }
    for (;;) {
        double lower_bound;
        tep_status = tep_generator_peek(&generator, &lower_bound);
        if (tep_status == TEP_DONE) {
            result->optimality_proven = true;
            result->stop_reason = DECODE_STOP_EXHAUSTED;
            status = DECODE_SUCCESS;
            if (trace) fputs("trace stop=exhausted (no remaining f_P)\n", trace);
            break;
        }
        if (tep_status != TEP_OK) goto tep_error;
        if (result->candidate_valid && lower_bound >= result->cost) {
            result->optimality_proven = true;
            result->stop_reason = DECODE_STOP_BOUND;
            status = DECODE_SUCCESS;
            if (trace) fprintf(trace, "trace stop=bound next_W_P=%.17g >= best=%.17g\n", lower_bound, result->cost);
            break;
        }
        if (config->candidate_limit && result->stats.candidates_evaluated >= config->candidate_limit) {
            decode_result_mark_truncated(result);
            status = DECODE_TRUNCATED;
            if (trace) fprintf(trace, "trace stop=budget evaluated=%" PRIu64 " limit=%" PRIu64
                               " next_W_P=%.17g < best=%.17g\n", result->stats.candidates_evaluated,
                               config->candidate_limit, lower_bound, result->cost);
            break;
        }
        double cost;
        tep_status = tep_generator_next(&generator, &partial, &cost);
        if (tep_status != TEP_OK) goto tep_error;
        double partial_cost = cost;
        ++result->stats.candidates_popped;
        memset(system.data, 0, m);
        if (k) memcpy(system.data + m, partial.data, k);
        if (m) {
            /* H_sys * (0,f_P)^T = P*f_P^T; views reuse S03 storage/scratch. */
            BinaryMatrix h_sys = {m, n, form->reduced.data};
            BitVector basis_syndrome = {m, scratch.data};
            if (binary_matrix_multiply_vector(&h_sys, &system, &basis_syndrome) != GF2_OK) {
                status = DECODE_INVALID_ARGUMENT;
                goto finish;
            }
            for (size_t i = 0; i < m; ++i) system.data[i] = transformed.data[i] ^ scratch.data[i];
        }
        /* W = W_P + W_I, with only nonnegative additions: W_P remains a valid
         * floating-point lower bound. This is the same flip objective W(f). */
        for (size_t i = 0; i < m; ++i) {
            if (system.data[i]) cost += weights[form->sys_to_orig[i]];
        }
        if (!isfinite(cost)) { status = DECODE_NUMERIC_ERROR; goto finish; }
        ++result->stats.candidates_evaluated;
        bool improved = !result->candidate_valid || cost < result->cost;
        if (improved) {
            if (binary_systematic_to_original(form, &system, &original) != GF2_OK) {
                status = DECODE_INVALID_ARGUMENT;
                goto finish;
            }
            for (size_t j = 0; j < n; ++j) result->estimate.data[j] = reference.data[j] ^ original.data[j];
            result->cost = cost;
            result->candidate_valid = true;
        }
        if (trace) {
            fprintf(trace, "trace candidate=%" PRIu64, result->stats.candidates_evaluated);
            trace_bits(trace, " f_P(sys)=", partial.data, k);
            fprintf(trace, " W_P=%.17g", partial_cost);
            trace_bits(trace, " f_I(sys)=", system.data, m);
            fprintf(trace, " W=%.17g best=%.17g improved=%d\n", cost, result->cost, (int)improved);
        }
    }
    goto finish;
tep_error:
    status = tep_status == TEP_ALLOCATION_FAILED ? DECODE_ALLOCATION_FAILED :
             tep_status == TEP_NUMERIC_ERROR ? DECODE_NUMERIC_ERROR : DECODE_INVALID_ARGUMENT;
finish:
    result->stats.nodes_expanded = generator.nodes_expanded;
    if (result->candidate_valid) {
        /* Verify the returned best estimate against the ORIGINAL H/target. */
        int valid = binary_matrix_multiply_vector(input->h, &result->estimate, &scratch) == GF2_OK;
        for (size_t i = 0; valid && i < rows; ++i) {
            uint8_t target = input->mode == DECODE_CLASSICAL ? 0 : input->s_target->data[i];
            valid = scratch.data[i] == target;
        }
        if (!valid) {
            status = DECODE_INVALID_CANDIDATE;
            result->candidate_valid = result->optimality_proven = false;
            result->stop_reason = DECODE_STOP_NONE;
        }
    }
    tep_generator_free(&generator);
    for (size_t i = 0; i < sizeof(buffers) / sizeof(*buffers); ++i) bit_vector_free(buffers[i]);
    free(weights);
    free(partial_weights);
    return result->status = status;
}

bool decoder_is_known(const char *name)
{
    static const char *const names[] = {"gcd", "lcgcd", "bp", "nms", "bp-lcgcd", "bp-osd"};
    if (!name) return false;
    for (size_t i = 0; i < sizeof(names) / sizeof(names[0]); ++i) {
        if (strcmp(name, names[i]) == 0) return true;
    }
    return false;
}

//译码--> gcd_decode
DecodeStatus decoder_decode(const DecodeConfig *config, const DecodeInput *input,
                            DecodeResult *result)
{
    if (!result) return DECODE_INVALID_ARGUMENT;
    decode_result_reset(result);
    if (!config || !decoder_is_known(config->decoder_name)) return result->status = DECODE_INVALID_ARGUMENT;
    if (strcmp(config->decoder_name, "gcd") == 0) return gcd_decode(config, input, result);
    DecodeStatus status = validate_input(input);
    if (status != DECODE_SUCCESS) return result->status = status;
    if (!result->estimate.data) {
        return result->status = DECODE_INVALID_ARGUMENT;
    }
    if (result->estimate.length != input->length) {
        return result->status = DECODE_DIMENSION_MISMATCH;
    }
    return result->status = DECODE_NOT_IMPLEMENTED;
}
