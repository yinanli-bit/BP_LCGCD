#include "bplcgcd/decoder.h"
#include "bplcgcd/tep.h"
#include "bplcgcd/slva.h"

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
    
    /* delta=0 reduces exactly to GCD, including TEP ties and budget counts. */
    bool lc = strcmp(config->decoder_name, "lcgcd") == 0 && config->delta != 0;
    size_t k = n - m, delta = lc ? config->delta : 0;
    if (delta > m) 
        return result->status = DECODE_DIMENSION_MISMATCH;
    if (lc && ((m && k > SIZE_MAX/m) || n > SIZE_MAX/sizeof(size_t)))
        return result->status = DECODE_ALLOCATION_FAILED;
    size_t stages = lc ? delta+k : k;
    BitVector reference = {0}, residual = {0}, transformed = {0}, system = {0};
    BitVector partial = {0}, original = {0}, scratch = {0}, local_rhs = {0};
    BitVector *buffers[] = {&reference, &residual, &transformed, &system, &partial, &original, &scratch, &local_rhs};
    const size_t lengths[] = {n, rows, rows, n, stages, n, rows, lc ? m : 0};
    TepGenerator generator = {0};

    SlvaGenerator list = {0};
    SyndromeTrellis trellis = {.length=stages};
    LocalConstraints local = {0};

    uint8_t *local_p = NULL;
    TepStatus tep_status = TEP_OK;
    SlvaStatus slva_status = SLVA_OK;
    DecodeStatus status = DECODE_ALLOCATION_FAILED;
    double *weights = malloc(n * sizeof(double));
    double *partial_weights = !lc && k ? malloc(k * sizeof(double)) : NULL;
    if (!weights || (!lc && k && !partial_weights)) goto finish;
    for (size_t i = 0; i < sizeof(buffers) / sizeof(*buffers); ++i) {
        if (lengths[i] && bit_vector_init(buffers[i], lengths[i]) != GF2_OK) goto finish;
    }
    status = decode_prepare(input, &reference, &residual, weights);
    if (status != DECODE_SUCCESS) goto finish;
    Gf2Status gf2_status;
    if (lc) {
        /* The two P blocks and right-hand sides are views into owned buffers. */
        local_p = m && k ? malloc(m*k) : NULL;
        local.basis_order = m ? malloc(m*sizeof(size_t)) : NULL;
        trellis.column_masks = stages ? malloc(stages*sizeof(size_t)) : NULL;
        trellis.original_columns = stages ? malloc(stages*sizeof(size_t)) : NULL;
        if ((m && k && !local_p) || (m && !local.basis_order) ||
            (stages && (!trellis.column_masks || !trellis.original_columns))) {
            status = DECODE_ALLOCATION_FAILED;
            goto finish;
        }
        local.p_a = (BinaryMatrix){delta,k,local_p};
        local.p_ac = (BinaryMatrix){m-delta,k,local_p ? local_p+delta*k : NULL};
        local.s_a = (BitVector){delta,local_rhs.data};
        local.s_ac = (BitVector){m-delta,local_rhs.data ? local_rhs.data+delta : NULL};
        gf2_status = local_constraints_prepare(form, weights, n, &residual, delta, &local, &transformed);
    } 
    else {
        gf2_status = binary_systematic_transform_syndrome(form, &residual, &transformed);
    }
    if (gf2_status != GF2_OK) {
        status = gf2_status == GF2_INCONSISTENT_SYNDROME ?
                 DECODE_INCONSISTENT_SYNDROME : DECODE_INVALID_ARGUMENT;
        goto finish;
    }
    if (lc) {
        TrellisStatus graph_status = syndrome_trellis_prepare(form, &local, weights, n,config->state_limit ? config->state_limit : 4096, &trellis);
        if (graph_status != TRELLIS_OK) {
            slva_status = graph_status == TRELLIS_RESOURCE_LIMIT ? SLVA_RESOURCE_LIMIT :
                          graph_status == TRELLIS_NUMERIC_ERROR ? SLVA_NUMERIC_ERROR : SLVA_INVALID_ARGUMENT;
            goto slva_error;
        }
        SlvaLimits limits = {.max_candidates=config->slva_max_candidates ? config->slva_max_candidates : 65536,.max_bytes=config->slva_max_bytes ? config->slva_max_bytes : 64*1024*1024};
        slva_status = slva_init(&list, &trellis, &limits);
        if (slva_status != SLVA_OK) goto slva_error;
    } else {
        for (size_t j = 0; j < k; ++j) partial_weights[j] = weights[form->sys_to_orig[m+j]];
        tep_status = tep_generator_init(&generator, partial_weights, k);
        if (tep_status != TEP_OK) goto tep_error;
    }
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
        if (lc) {
            fprintf(trace, "trace local delta=%zu A(original)=", delta);
            for (size_t i = 0; i < delta; ++i) fprintf(trace, "%s%zu", i ? "," : "", local.basis_order[i]);
            trace_bits(trace, " s_A=", local.s_a.data, delta);
            fputc('\n', trace);
        }
    }
    double partial_cost = 0;
    for (;;) {
        /* All unseen local paths cost >= the last returned one. In particular
         * delta=m needs only Viterbi-1, with no SLVA expansion for a proof. */
        if (lc && result->candidate_valid && partial_cost >= result->cost) {
            result->optimality_proven = true;
            result->stop_reason = DECODE_STOP_BOUND;
            status = DECODE_SUCCESS;
            if (trace) fprintf(trace, "trace stop=bound last_W_local=%.17g >= best=%.17g\n", partial_cost, result->cost);
            break;
        }
        double lower_bound;
        if (lc) {
            slva_status = slva_peek(&list, &lower_bound);
            if (slva_status != SLVA_OK && slva_status != SLVA_DONE) goto slva_error;
        } else {
            tep_status = tep_generator_peek(&generator, &lower_bound);
            if (tep_status != TEP_OK && tep_status != TEP_DONE) goto tep_error;
        }
        if (lc ? slva_status == SLVA_DONE : tep_status == TEP_DONE) {
            result->optimality_proven = true;
            result->stop_reason = DECODE_STOP_EXHAUSTED;
            status = DECODE_SUCCESS;
            if (trace) fprintf(trace, "trace stop=exhausted (no remaining %s)\n", lc ? "local paths" : "f_P");
            break;
        }
        if (result->candidate_valid && lower_bound >= result->cost) {
            result->optimality_proven = true;
            result->stop_reason = DECODE_STOP_BOUND;
            status = DECODE_SUCCESS;
            if (trace) fprintf(trace, "trace stop=bound %s=%.17g >= best=%.17g\n",
                               lc ? "next_W_local" : "next_W_P", lower_bound, result->cost);
            break;
        }
        if (config->candidate_limit && result->stats.candidates_evaluated >= config->candidate_limit) {
            decode_result_mark_truncated(result);
            status = DECODE_TRUNCATED;
            if (trace) fprintf(trace, "trace stop=budget evaluated=%" PRIu64 " limit=%" PRIu64
                               " %s=%.17g < best=%.17g\n", result->stats.candidates_evaluated,
                               config->candidate_limit, lc ? "next_W_local" : "next_W_P", lower_bound, result->cost);
            break;
        }
        double cost;
        if (lc) {
            slva_status = slva_next(&list, &partial, &cost);
            if (slva_status != SLVA_OK) goto slva_error;
        } else {
            tep_status = tep_generator_next(&generator, &partial, &cost);
            if (tep_status != TEP_OK) goto tep_error;
        }
        partial_cost = cost;
        ++result->stats.candidates_popped;

        // 补全完整的system，lcgcd分三步（k,delta,m-delta）,gcd分两步（k,m）
        memset(system.data, 0, m);
        if (k) memcpy(system.data + m, partial.data + (lc ? delta : 0), k); // copy the last k bits
        if (lc) {
            for (size_t i = 0; i < delta; ++i)
                system.data[form->orig_to_sys[local.basis_order[i]]] = partial.data[i]; //按照basis_order的顺序，把前delta个bit放到system中
            BitVector free_bits = {k,system.data+m};
            BitVector product = {m-delta,scratch.data}; //需要补全的bit
            //P_Ac*f_Pᵀ
            if (m > delta && k && binary_matrix_multiply_vector(&local.p_ac, &free_bits, &product) != GF2_OK) {
                status = DECODE_INVALID_ARGUMENT;
                goto finish;
            }
            for (size_t i = 0; i < m-delta; ++i)
                //f_I[Ac] = s_Ac XOR P_Ac*f_Pᵀ 补齐候选
                system.data[form->orig_to_sys[local.basis_order[delta+i]]] = local.s_ac.data[i] ^ (k ? scratch.data[i] : 0); 
        } else if (m) {
            /* H_sys * (0,f_P)^T = P*f_P^T; views reuse S03 storage/scratch. */
            BinaryMatrix h_sys = {m, n, form->reduced.data};
            BitVector basis_syndrome = {m, scratch.data};
            if (binary_matrix_multiply_vector(&h_sys, &system, &basis_syndrome) != GF2_OK) {
                status = DECODE_INVALID_ARGUMENT;
                goto finish;
            }
            for (size_t i = 0; i < m; ++i) 
                //e_I=s XOR e_P P^T 补齐候选
                system.data[i] = transformed.data[i] ^ scratch.data[i];
        }
        //计算完整代价
        /* GCD: W_P+W_I; LC-GCD: W_local+W_Ac. Append nonnegative terms to
         * the enumerator's score, preserving its floating-point lower bound. */
        for (size_t i = lc ? delta : 0; i < m; ++i) {
            size_t original_column = lc ? local.basis_order[i] : form->sys_to_orig[i];
            if (system.data[form->orig_to_sys[original_column]]) 
                cost += weights[original_column];
        }
        if (!isfinite(cost)) { status = DECODE_NUMERIC_ERROR; goto finish; }
        if (lc) {
            if (binary_systematic_to_original(form, &system, &original) != GF2_OK ||binary_matrix_multiply_vector(input->h, &original, &scratch) != GF2_OK) {
                status = DECODE_INVALID_ARGUMENT;
                goto finish;
            }
            for (size_t i = 0; i < rows; ++i) if (scratch.data[i] != residual.data[i]) {
                status = DECODE_INVALID_CANDIDATE;
                goto finish;
            }
        }
        ++result->stats.candidates_evaluated;

        // Update the best candidate.
        bool improved = !result->candidate_valid || cost < result->cost;
        if (improved) {
            //下面这一个判断条件，lcgcd跳过，gcd则计算binary_systematic_to_original
            if (!lc && binary_systematic_to_original(form, &system, &original) != GF2_OK) {
                status = DECODE_INVALID_ARGUMENT;
                goto finish;
            }
            for (size_t j = 0; j < n; ++j) 
                result->estimate.data[j] = reference.data[j] ^ original.data[j]; //得到码字(classic)或错误模式(the other two)
            result->cost = cost; //更新
            result->candidate_valid = true;
        }
        if (trace) {
            fprintf(trace, "trace candidate=%" PRIu64, result->stats.candidates_evaluated);
            trace_bits(trace, lc ? " f_A/P(stage)=" : " f_P(sys)=", partial.data, stages);
            fprintf(trace, " %s=%.17g", lc ? "W_local" : "W_P", partial_cost);
            trace_bits(trace, " f_I(sys)=", system.data, m);
            fprintf(trace, " W=%.17g best=%.17g improved=%d\n", cost, result->cost, (int)improved);
        }
    }
    goto finish;
slva_error:
    if (slva_status == SLVA_RESOURCE_LIMIT || slva_status == SLVA_PATH_LIMIT) {
        decode_result_mark_truncated(result);
        result->stop_reason = DECODE_STOP_RESOURCE;
        status = DECODE_TRUNCATED;
        if (config->trace) fputs("trace stop=resource (Trellis/SLVA cap; no optimality proof)\n", config->trace);
    } else status = slva_status == SLVA_ALLOCATION_FAILED ? DECODE_ALLOCATION_FAILED :
                    slva_status == SLVA_NUMERIC_ERROR ? DECODE_NUMERIC_ERROR : DECODE_INVALID_ARGUMENT;
    goto finish;
tep_error:
    status = tep_status == TEP_ALLOCATION_FAILED ? DECODE_ALLOCATION_FAILED :
             tep_status == TEP_NUMERIC_ERROR ? DECODE_NUMERIC_ERROR : DECODE_INVALID_ARGUMENT;
finish:
    result->stats.nodes_expanded = lc ? list.nodes_expanded : generator.nodes_expanded;
    // next_sequence = 根代表 + 已成功提交的 SLVA 子组代表。
    result->stats.partial_patterns_generated = lc ? list.next_sequence : generator.partial_patterns_generated;
    if (lc && list.cells) {
        for (size_t stage=0; stage<trellis.length; ++stage) {
            for (size_t state=0; state<trellis.state_count; ++state) {
                if (isfinite(list.cells[stage*trellis.state_count+state].metric))
                    ++result->stats.trellis_states_visited;
            }
        }
    }
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
    slva_free(&list);
    free(local_p);
    free(local.basis_order);
    free(trellis.column_masks);
    free(trellis.original_columns);
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
    if (strcmp(config->decoder_name, "gcd") == 0 || strcmp(config->decoder_name, "lcgcd") == 0)
        return gcd_decode(config, input, result);
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
