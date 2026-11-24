#ifndef BPLCGCD_VITERBI_H
#define BPLCGCD_VITERBI_H

#include "bplcgcd/trellis.h"

/* Caller-owned full prefix table: cells[layer*state_count + state].
 * Unreachable cells have metric=INFINITY, previous_state=SIZE_MAX.
 * The bit label is essential for parallel edges with identical endpoints. */
typedef struct {
    double metric; //到达当前状态的最小累计代价
    size_t previous_state; //最优路径来自上一层的哪个状态
    uint8_t bit; //到达当前状态时选择的是 0 还是 1
} ViterbiCell;


TrellisStatus viterbi_decode(const SyndromeTrellis *trellis, ViterbiCell *cells,
                            size_t cell_capacity, BitVector *path, double *cost);

#endif
