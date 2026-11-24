#include "bplcgcd/trellis.h"

#include <limits.h>
#include <math.h>

void print_trellis(SyndromeTrellis *trellis) {
    printf("trellis->length = %zu\n", trellis->length);
    //printf("trellis->state_limit = %zu\n", trellis->state_limit);
    printf("trellis->flip_weights = [");
    for (size_t i = 0; i < trellis->length; ++i) {
        printf("%f", trellis->flip_weights[i]);
        if (i < trellis->length - 1) {
            printf(", ");
        }
    }
    printf("]\n");
    printf("trellis->column_masks = [");
    for (size_t i = 0; i < trellis->length; ++i) {
        printf("%zu", trellis->column_masks[i]);
        if (i < trellis->length - 1) {
            printf(", ");
        }
    }
    printf("]\n");
    printf("trellis->original_columns = [");
    for (size_t i = 0; i < trellis->length; ++i) {
        printf("%zu", trellis->original_columns[i]);
        if (i < trellis->length - 1) {
            printf(", ");
        }
    }
    printf("]\n");
}

//准备图
TrellisStatus syndrome_trellis_prepare(const BinarySystematic *form, const LocalConstraints *local,
                                      const double *flip_weights, size_t length, size_t state_limit,
                                      SyndromeTrellis *result)
{
    if (!form || !local || !flip_weights || !result || !form->sys_to_orig || !state_limit) {
        return TRELLIS_INVALID_ARGUMENT;
    }
    size_t n = form->reduced.cols, m = form->rank, delta = local->p_a.rows;
    if (m > n || delta > m || length != n) 
        return TRELLIS_DIMENSION_MISMATCH;
    
    size_t k = n - m;
    size_t stages = delta + k; //δ+k

    if (local->p_a.cols != k || local->s_a.length != delta || result->length != stages) {
        return TRELLIS_DIMENSION_MISMATCH;
    }
    if (delta >= sizeof(size_t) * CHAR_BIT) 
        return TRELLIS_RESOURCE_LIMIT;
    
    size_t states = (size_t)1 << delta; //2^δ

    if (states > state_limit) 
        return TRELLIS_RESOURCE_LIMIT;
    if ((stages && (!result->column_masks || !result->original_columns)) ||(delta && (!local->basis_order || !local->s_a.data)) ||(delta && k && !local->p_a.data)) 
        return TRELLIS_INVALID_ARGUMENT;

    size_t target = 0;//目标 syndrome 的整数表示
    for (size_t i = 0; i < delta; ++i) 
        target |= (size_t)local->s_a.data[i] << i;
    double total = 0;
    //column_masks 记录每个 Trellis 比特对局部 syndrome 的影响；
    //original_columns 记录这个比特实际对应原始校验矩阵中的哪一列。
    for (size_t stage = 0; stage < stages; ++stage) {
        size_t mask = 0, original;
        if (stage < delta) {
            mask = (size_t)1 << stage;
            original = local->basis_order[stage];
        } else {
            size_t j = stage - delta;
            for (size_t i = 0; i < delta; ++i) 
                mask |= (size_t)local->p_a.data[i*k+j] << i;
            original = form->sys_to_orig[m+j];
        }
        double weight = flip_weights[original];
        if (!isfinite(weight) || weight < 0) return TRELLIS_INVALID_ARGUMENT;
        total += weight;
        if (!isfinite(total)) return TRELLIS_NUMERIC_ERROR;
        result->column_masks[stage] = mask;
        result->original_columns[stage] = original;
    }
    result->delta = delta;
    result->state_count = states;
    result->target = target;
    result->flip_weights = flip_weights;

    //print_trellis(result);

    return TRELLIS_OK;
}

TrellisStatus syndrome_trellis_branch(const SyndromeTrellis *trellis, size_t stage, size_t state,
                                     uint8_t bit, size_t *next_state, double *branch_cost)
{
    if (!trellis || !next_state || !branch_cost || bit > 1) 
        return TRELLIS_INVALID_ARGUMENT;
    if (stage >= trellis->length || state >= trellis->state_count) 
        return TRELLIS_DIMENSION_MISMATCH;
    *next_state = state ^ (bit ? trellis->column_masks[stage] : 0);
    *branch_cost = bit ? trellis->flip_weights[trellis->original_columns[stage]] : 0;

     //printf("stage: %zu, state: %zu, bit: %u, next_state: %zu, branch_cost: %f\n", stage, state, bit, *next_state, *branch_cost);
    return TRELLIS_OK;
}
