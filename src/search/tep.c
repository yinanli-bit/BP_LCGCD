#include "bplcgcd/tep.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct TepNode {
    double cost;
    size_t hw;
    size_t *indices; /* Selected sorted positions, increasing and zero-based. */
} TepNode;

static int before(const TepNode *a, const TepNode *b)
{
    if (a->cost != b->cost) return a->cost < b->cost;
    size_t i = a->hw, j = b->hw;
    while (i && j) {
        size_t a_bit = a->indices[--i], b_bit = b->indices[--j];
        if (a_bit != b_bit) return a_bit < b_bit;
    }
    return i < j;
}

static TepStatus reserve(TepGenerator *g, size_t needed)
{
    if (needed <= g->capacity) return TEP_OK;
    size_t capacity = g->capacity ? 2 * g->capacity : 8;
    if (capacity > SIZE_MAX / sizeof(TepNode)) return TEP_ALLOCATION_FAILED;
    TepNode *nodes = realloc(g->nodes, capacity * sizeof(*nodes));
    if (!nodes) return TEP_ALLOCATION_FAILED;
    g->nodes = nodes;
    g->capacity = capacity;
    return TEP_OK;
}

//插入新候选，新候选上浮，比cost
static void push(TepGenerator *g, TepNode node)
{
    ++g->partial_patterns_generated;
    size_t i = g->size++; //先赋值 后自增
    while (i && before(&node, &g->nodes[(i - 1) / 2])) {
        g->nodes[i] = g->nodes[(i - 1) / 2];
        i = (i - 1) / 2;
    }
    g->nodes[i] = node;
}

// 取出最小cost候选，候选下沉，比cost
static TepNode pop(TepGenerator *g)
{
    TepNode root = g->nodes[0], tail = g->nodes[--g->size];
    size_t i = 0;
    while (2 * i + 1 < g->size) {
        size_t child = 2 * i + 1;
        if (child + 1 < g->size && before(&g->nodes[child + 1], &g->nodes[child])) ++child;
        if (!before(&g->nodes[child], &tail)) break;
        g->nodes[i] = g->nodes[child];
        i = child;
    }
    if (g->size) g->nodes[i] = tail;
    return root;
}

void tep_generator_free(TepGenerator *g)
{
    if (g) {
        for (size_t i = 0; i < g->size; ++i) free(g->nodes[i].indices);
        free(g->nodes);
        free(g->order);
        *g = (TepGenerator){0};
    }
}

//g 从{0}初始化；校验权重并建立内部排序和空根节点；失败保持 g 为空
TepStatus tep_generator_init(TepGenerator *g, const double *weights, size_t length)
{
    if (!g || (length && !weights) || length > SIZE_MAX / sizeof(size_t)) return TEP_INVALID_ARGUMENT;
    for (size_t j = 0; j < length; ++j) {
        if (!isfinite(weights[j]) || weights[j] < 0) return TEP_INVALID_ARGUMENT;
    }
    TepGenerator built = {.length = length, .flip_weights = weights}; //结构体全部初始化为0，length和flip_weights初始化为输入参数
    if (length) {
        built.order = malloc(length * sizeof(*built.order));
        if (!built.order) return TEP_ALLOCATION_FAILED;
    }
    /* Stable insertion sort: equal weights retain increasing input index. */
    // 按权重递增排序，排序后的位置映射到原始比特位置
    double total = 0;
    for (size_t j = 0; j < length; ++j) {
        size_t i = j;
        while (i && weights[built.order[i - 1]] > weights[j]) {
            built.order[i] = built.order[i - 1];
            --i;
        }
        built.order[i] = j;
    }
    printf("order: ");
    for (size_t j = 0; j < length; ++j) printf("%zu ", built.order[j]);
    printf("\n");

    for (size_t j = 0; j < length; ++j) total += weights[built.order[j]];
    if (!isfinite(total)) {
        tep_generator_free(&built);
        return TEP_NUMERIC_ERROR;
    }
    TepStatus status = reserve(&built, 1);
    if (status != TEP_OK) {
        tep_generator_free(&built);
        return status;
    }
    push(&built, (TepNode){0}); /* Empty pattern, hw=0, no selected positions. */
    *g = built;
    return TEP_OK;
}
void print_heap(const TepGenerator *g) {
    if (!g) {
        printf("Heap is NULL\n");
        return;
    }
    printf("Heap size: %zu\n", g->size);
    for (size_t i = 0; i < g->size; ++i) {
        TepNode node = g->nodes[i];
        printf("Node %zu: cost=%.6f, hw=%zu, indices=[", i, node.cost, node.hw);
        for (size_t j = 0; j < node.hw; ++j) {
            printf("%zu", node.indices[j]);
            if (j + 1 < node.hw) printf(", ");
        }
        printf("]\n");
    }
}
//每个TEP_OK恰好一次出堆、一个候选，包括空模式，一次取一个最优出来
TepStatus tep_generator_next(TepGenerator *g, BitVector *pattern, double *cost)
{
    if (!g || !pattern || !cost || (g->length && !pattern->data)) return TEP_INVALID_ARGUMENT;
    if (pattern->length != g->length) return TEP_DIMENSION_MISMATCH;
    if (!g->size) return TEP_DONE;
    TepNode parent = g->nodes[0], children[2] = {0};
    size_t next = parent.hw ? parent.indices[parent.hw - 1] + 1 : 0;
    size_t count = next < g->length ? (parent.hw ? 2 : 1) : 0;
    /* Prepare both children and capacity before popping, so allocation errors
     * leave the iterator at the same pattern. No full subset list is built. */
    for (size_t c = 0; c < count; ++c) {
        //The empty parent has only the special right child {0}. 
        int right = !parent.hw || c == 1;
        children[c].hw = parent.hw + right;
        children[c].indices = malloc(children[c].hw * sizeof(size_t));
        if (!children[c].indices) {
            if (c) free(children[0].indices);
            return TEP_ALLOCATION_FAILED;
        }
        if (parent.hw) memcpy(children[c].indices, parent.indices, parent.hw * sizeof(size_t));
        children[c].indices[children[c].hw - 1] = next;
        /* Sum the selected sorted weights directly, avoiding subtraction-based
         * cost updates and preserving the parent <= child ordering in double. */
        for (size_t j = 0; j < children[c].hw; ++j) {
            children[c].cost += g->flip_weights[g->order[children[c].indices[j]]];
        }
    }
    TepStatus status = reserve(g, g->size - 1 + count);
    if (status != TEP_OK) {
        for (size_t c = 0; c < count; ++c) free(children[c].indices);
        return status;
    }

    
    printf("当前二叉堆的内部状态：\n");
    print_heap(g);
    

    parent = pop(g); //取出最优候选
    for (size_t c = 0; c < count; ++c) push(g, children[c]);
    if (count) ++g->nodes_expanded;
    if (g->length) memset(pattern->data, 0, g->length);
    for (size_t j = 0; j < parent.hw; ++j) pattern->data[g->order[parent.indices[j]]] = 1;
    *cost = parent.cost;
    free(parent.indices);
    return TEP_OK;
}

// 只读查询剩余最小部分代价，不分配、出堆或扩展
TepStatus tep_generator_peek(const TepGenerator *g, double *cost)
{
    if (!g || !cost) return TEP_INVALID_ARGUMENT;
    if (!g->size) return TEP_DONE;
    *cost = g->nodes[0].cost;
    return TEP_OK;
}

/*FPT 树：规定一个候选如何生成后续候选。
  最小堆：保存已经生成但尚未输出的候选，并保证下一次 pop() 取出其中代价最小的。*/
