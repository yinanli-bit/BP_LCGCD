#ifndef BPLCGCD_SLVA_H
#define BPLCGCD_SLVA_H

#include "bplcgcd/viterbi.h"

typedef enum {
    SLVA_OK = 0,
    SLVA_DONE,
    SLVA_PATH_LIMIT,
    SLVA_RESOURCE_LIMIT,
    SLVA_ALLOCATION_FAILED,
    SLVA_INVALID_ARGUMENT,
    SLVA_DIMENSION_MISMATCH,
    SLVA_NUMERIC_ERROR
} SlvaStatus;

typedef struct {
    size_t max_candidates; /* Live heap/pending/temporary paths; must be >0. */
    size_t max_bytes;      /* Owned allocation payload, including expansion; >0. */
    uint64_t path_limit;   /* Returned paths, 0=unlimited (resource caps still apply). */
} SlvaLimits;

/* Internal iterator state, read-only to callers. A candidate represents all
 * paths with bits[free_length..L) fixed; its prefix is a Viterbi survivor. */
typedef struct {
    uint8_t *bits; //当前候选组的最轻代表路径，按 Trellis 阶段顺序保存
    size_t free_length; //此分组中允许变化的前缀长度
    double cost; //该代表路径的局部代价
    uint64_t sequence; //候选创建顺序
} SlvaCandidate;

typedef struct {
    const SyndromeTrellis *trellis; //指向 S11 构造的 Trellis，包含状态转移信息、各阶段权重和目标 syndrome
    ViterbiCell *cells; //指向 S12 构造的 Viterbi 表，包含每个阶段每个状态的最优前驱、最优比特和最优累计代价
    SlvaCandidate *heap; //最小堆，存放候选组及其最轻代表路径
    size_t size, capacity; //当前堆中实际有多少个候选组，以及堆的容量
    SlvaCandidate pending; //已经输出、但剩余路径还没有拆分成子组的候选
    SlvaLimits limits; //用户设置的候选数量、内存和输出路径数量限制
    size_t bytes_used; //当前 SLVA 模块管理的内存占用量
    uint64_t paths_returned, next_sequence; //已经通过 slva_next() 输出的路径数量，为新候选分配的顺序号
} SlvaGenerator;


SlvaStatus slva_init(SlvaGenerator *generator, const SyndromeTrellis *trellis, const SlvaLimits *limits);

SlvaStatus slva_next(SlvaGenerator *generator, BitVector *path, double *cost);

/* Does not emit a path; may lazily expand the previous output and allocate.
 * Repeated successful peeks return the same next cost. */
SlvaStatus slva_peek(SlvaGenerator *generator, double *cost);
void slva_free(SlvaGenerator *generator);

#endif
