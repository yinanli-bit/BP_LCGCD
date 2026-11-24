#include "bplcgcd/viterbi.h"

#include <math.h>

void print_viterbi_cells(const ViterbiCell *cells, size_t length, size_t state_count) {
    for (size_t layer = 0; layer <= length; ++layer) {
        printf("Layer %zu:\n", layer);
        for (size_t state = 0; state < state_count; ++state) {
            const ViterbiCell *cell = &cells[layer * state_count + state];
            printf("  State %zu: metric=%f, previous_state=%zu, bit=%u\n",
                   state, cell->metric, cell->previous_state, cell->bit);
        }
    }
}
//找到最小代价路径，返回路径和代价：找到可能的前驱-->前缀代价 + 分支代价-->保留更优的路径
TrellisStatus viterbi_decode(const SyndromeTrellis *trellis, ViterbiCell *cells,
                            size_t cell_capacity, BitVector *path, double *cost)
{
    if (!trellis || !cells || !path || !cost || !trellis->state_count ||trellis->target >= trellis->state_count) 
        return TRELLIS_INVALID_ARGUMENT;
    
    size_t length = trellis->length, states = trellis->state_count;

    if (path->length != length) 
        return TRELLIS_DIMENSION_MISMATCH;
    if (length && (!path->data || !trellis->column_masks || !trellis->original_columns ||!trellis->flip_weights)) 
        return TRELLIS_INVALID_ARGUMENT;
    if (length == SIZE_MAX || states > SIZE_MAX / (length+1)) 
        return TRELLIS_RESOURCE_LIMIT;
    
    size_t required = (length+1)*states;

    if (required > SIZE_MAX / sizeof(*cells) || required > cell_capacity) 
        return TRELLIS_RESOURCE_LIMIT;

    //cells初始化
    for (size_t i = 0; i < required; ++i) {
        cells[i] = (ViterbiCell){.metric=INFINITY, .previous_state=SIZE_MAX};
    }
    cells[0].metric = 0;
    /*
    for(size_t stage = 0; stage < length; ++stage){
        printf("column mask: %zu\n", trellis->column_masks[stage]);
        printf("original column: %zu\n", trellis->original_columns[stage]);
        printf("flip weight: %f\n", trellis->flip_weights[stage]);
        
    }*/
    for (size_t stage = 0; stage < length; ++stage) {
        for (size_t state = 0; state < states; ++state) {
            ViterbiCell *best = &cells[(stage+1)*states+state];//指向下一层，best相当于cell的数组cells[],state从0变成1时，best指向的是另一个格子，并不会发生覆盖
            for (uint8_t bit = 0; bit < 2; ++bit) {
                /* XOR is its own inverse: the two incoming labeled edges. */
                size_t previous = state ^ (bit ? trellis->column_masks[stage] : 0);
                double prefix = cells[stage*states+previous].metric; //对应于前一个格子的best->metric
                //printf("prefix: %f\n", prefix);
                if (!isfinite(prefix)) continue;
                size_t next;
                double edge;
                TrellisStatus status = syndrome_trellis_branch(trellis, stage, previous, bit, &next, &edge);//返回后继和本阶段非负增量代价
                if (status != TRELLIS_OK) 
                    return status;
                
                double candidate = prefix + edge;

                if (!isfinite(candidate)) 
                    return TRELLIS_NUMERIC_ERROR;
                if (candidate < best->metric) {
                    //每个格子保存的是从初始状态出发的全局最优前缀
                    *best = (ViterbiCell){.metric=candidate, .previous_state=previous, .bit=bit};
                    //printf("===第%zu层找到的最优路径=== \n", stage+1);
                    //printf("stage: %zu, state: %zu, previous_state: %zu, bit: %u, metric: %f\n", stage, state, previous, bit, candidate);
                    //print_viterbi_cells(best, stage+1, states);
                }
            }
        }
        //printf("===第%zu层找到的最优路径=== \n", stage+1);
    }
    double optimum = cells[length*states+trellis->target].metric;
    if (!isfinite(optimum)) return TRELLIS_NO_PATH;
    size_t state = trellis->target;
    //回溯路径
    for (size_t layer = length; layer > 0; --layer) {
        const ViterbiCell *cell = &cells[layer*states+state];
        path->data[layer-1] = cell->bit;
        state = cell->previous_state;
        //printf("回溯第%zu层的最优路径：bit=%u, previous_state=%zu\n", layer, cell->bit, cell->previous_state);
    }
    printf("最优路径：\n");
    for (size_t i = 0; i < length; ++i) {
        printf("%u ", path->data[i]);
    }
    printf("\n");
    *cost = optimum;
    return TRELLIS_OK;
}
