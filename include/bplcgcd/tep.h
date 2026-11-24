#ifndef BPLCGCD_TEP_H
#define BPLCGCD_TEP_H

#include "bplcgcd/gf2.h"

typedef enum {
    TEP_OK = 0, // 成功，返回一个候选模式，对应一次出堆、一个候选
    TEP_DONE, // 搜索空间耗尽，返回空模式，表示没有更多候选
    TEP_INVALID_ARGUMENT,
    TEP_DIMENSION_MISMATCH,
    TEP_ALLOCATION_FAILED,
    TEP_NUMERIC_ERROR
} TepStatus;

// 二叉最小堆保存待扩展节点
typedef struct {
    size_t length; //错误模式的比特长度
    const double *flip_weights;
    size_t *order; //按权重递增排序后的原始比特位置映射
    struct TepNode *nodes;
    size_t size;
    size_t capacity;
    uint64_t nodes_expanded; //出堆后生成至少一个孩子的父节点数
} TepGenerator;


TepStatus tep_generator_init(TepGenerator *generator,const double *flip_weights, size_t length);

TepStatus tep_generator_next(TepGenerator *generator, BitVector *pattern,double *cost);
/* 只读查询堆顶代价；耗尽时返回 TEP_DONE，且不修改 cost。 */
TepStatus tep_generator_peek(const TepGenerator *generator, double *cost);
void tep_generator_free(TepGenerator *generator);

#endif
