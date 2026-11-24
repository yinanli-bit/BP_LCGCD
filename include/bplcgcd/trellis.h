#ifndef BPLCGCD_TRELLIS_H
#define BPLCGCD_TRELLIS_H

#include "bplcgcd/local_constraints.h"

typedef enum {
    TRELLIS_OK = 0,
    TRELLIS_INVALID_ARGUMENT,
    TRELLIS_DIMENSION_MISMATCH,
    TRELLIS_RESOURCE_LIMIT,
    TRELLIS_NUMERIC_ERROR,
    TRELLIS_NO_PATH
} TrellisStatus;

/*Trellis定义，列掩码 h_t 将矩阵第 t 列按相同位序打包，状态 q_t 是前 t 列已经累计出的 delta 位 syndrome
定义 `L=delta+K`、`M_local=[I_delta|P_A]`，第 t 个局部比特记 b_t。位序固定为 `(f_I[A],f_P)`*/
typedef struct {
    size_t length, delta, state_count, target;
    size_t *column_masks, *original_columns;
    const double *flip_weights;
} SyndromeTrellis;

//复用配对的 S03/S10 输出；weights 是原始坐标翻转权重，图使用期间保持不变
TrellisStatus syndrome_trellis_prepare(const BinarySystematic *form, const LocalConstraints *local,
                                      const double *flip_weights, size_t length, size_t state_limit,
                                      SyndromeTrellis *result);

/*  state XOR bit*column_masks[stage], at cost bit*weight(original_columns[stage]).
 * Both bit labels remain distinct even when they connect the same states. */
//返回后继和非负增量代价；不搜索路径
TrellisStatus syndrome_trellis_branch(const SyndromeTrellis *trellis, size_t stage, size_t state,
                                     uint8_t bit, size_t *next_state, double *branch_cost);

#endif
