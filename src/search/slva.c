#include "bplcgcd/slva.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <inttypes.h>

//slva_init() → slva_next() → alternative_count() → expand_pending() → push() / pop()


static SlvaStatus from_trellis(TrellisStatus status)
{
    switch (status) {
    case TRELLIS_OK: return SLVA_OK;
    case TRELLIS_DIMENSION_MISMATCH: return SLVA_DIMENSION_MISMATCH;
    case TRELLIS_RESOURCE_LIMIT: return SLVA_RESOURCE_LIMIT;
    case TRELLIS_NUMERIC_ERROR: return SLVA_NUMERIC_ERROR;
    default: return SLVA_INVALID_ARGUMENT;
    }
}

static void *allocate(SlvaGenerator *g, size_t bytes, SlvaStatus *status)
{
    if (bytes > g->limits.max_bytes-g->bytes_used) {
        *status = SLVA_RESOURCE_LIMIT;
        return NULL;
    }
    void *memory = malloc(bytes);
    if (!memory) { *status = SLVA_ALLOCATION_FAILED; return NULL; }
    g->bytes_used += bytes;
    *status = SLVA_OK;
    return memory;
}

static void release_path(SlvaGenerator *g, uint8_t *bits)
{
    if (bits) { free(bits); g->bytes_used -= g->trellis->length; }
}

/* Reserve room for heap nodes AND the path arrays about to be constructed. */
static SlvaStatus reserve(SlvaGenerator *g, size_t needed, size_t path_bytes)
{
    size_t available = g->limits.max_bytes-g->bytes_used;
    if (path_bytes > available) 
        return SLVA_RESOURCE_LIMIT;
    if (needed <= g->capacity) 
        return SLVA_OK;
    size_t capacity = g->capacity ? g->capacity : 1;
    if (capacity < needed) {
        capacity = capacity <= g->limits.max_candidates/2 ? 2*capacity : g->limits.max_candidates;
        if (capacity < needed) capacity = needed;
    }
    size_t budget = available-path_bytes;
    if (capacity > SIZE_MAX/sizeof(*g->heap) || (capacity-g->capacity)*sizeof(*g->heap) > budget) 
        capacity=needed;
    if (capacity > SIZE_MAX/sizeof(*g->heap) || (capacity-g->capacity)*sizeof(*g->heap) > budget) 
        return SLVA_RESOURCE_LIMIT;
    SlvaCandidate *heap = realloc(g->heap, capacity*sizeof(*heap));
    if (!heap) return SLVA_ALLOCATION_FAILED;
    g->bytes_used += (capacity-g->capacity)*sizeof(*heap);
    g->heap = heap;
    g->capacity = capacity;
    return SLVA_OK;
}

static int before(const SlvaCandidate *a, const SlvaCandidate *b)
{
    return a->cost < b->cost || (a->cost == b->cost && a->sequence < b->sequence);
}

//插入新候选，新候选上浮，比cost
static void push(SlvaGenerator *g, SlvaCandidate node)
{
    size_t i = g->size++;
    while (i && before(&node, &g->heap[(i-1)/2])) {
        g->heap[i] = g->heap[(i-1)/2];
        i = (i-1)/2;
    }
    g->heap[i] = node;
}

static SlvaCandidate pop(SlvaGenerator *g)
{
    SlvaCandidate root = g->heap[0], tail = g->heap[--g->size];
    size_t i = 0;
    while (i < g->size/2) {
        size_t child = 2*i+1;
        if (child+1 < g->size && before(&g->heap[child+1], &g->heap[child])) ++child;
        if (!before(&g->heap[child], &tail)) break;
        g->heap[i] = g->heap[child];
        i = child;
    }
    if (g->size) g->heap[i] = tail;
    return root;
}

void slva_free(SlvaGenerator *g)
{
    if (!g) return;
    for (size_t i = 0; i < g->size; ++i) free(g->heap[i].bits);
    free(g->pending.bits);
    free(g->heap);
    free(g->cells);
    *g = (SlvaGenerator){0};
}

static void print_slva_heap(SlvaCandidate *heap, size_t size){
    for (size_t i = 0; i < size; ++i) {
        printf("  %zu: %.17g\n", i, heap[i].cost);
    }
}
void print_slva_generator(const SlvaGenerator *g) {
    if (!g) {
        printf("SLVA Generator is NULL\n");
        return;
    }
    printf("SLVA Generator:\n");
    printf("  Size: %zu\n", g->size);
    printf("  Capacity: %zu\n", g->capacity);
    printf("  Bytes Used: %zu\n", g->bytes_used);
    printf("  Paths Returned: %" PRIu64 "\n", g->paths_returned);
    printf("  Next Sequence: %" PRIu64 "\n", g->next_sequence);
    printf("  Pending Candidate:\n");
    if (g->pending.bits) {
        printf("    Free Length: %zu\n", g->pending.free_length);
        printf("    Cost: %f\n", g->pending.cost);
        printf("    Sequence: %" PRIu64 "\n", g->pending.sequence);
        printf("    Bits: ");
        for (size_t j = 0; j < g->trellis->length; ++j) {
            printf("%u ", g->pending.bits[j]);
        }
        printf("\n");
    } else {
        printf("    No pending candidate.\n");
    }
}

//slva搜索器初始化，只调用一次完整 Viterbi，并保留全部前缀表供后续使用
SlvaStatus slva_init(SlvaGenerator *generator, const SyndromeTrellis *trellis, const SlvaLimits *limits)
{
    if (!generator || !trellis || !limits || !limits->max_candidates || !limits->max_bytes ||!trellis->state_count || trellis->target >= trellis->state_count) 
        return SLVA_INVALID_ARGUMENT;
    
    size_t length = trellis->length, states = trellis->state_count;

    if (length == SIZE_MAX || states > SIZE_MAX/(length+1)) 
        return SLVA_RESOURCE_LIMIT;

    size_t count = (length+1)*states;

    if (count > SIZE_MAX/sizeof(ViterbiCell)) 
        return SLVA_RESOURCE_LIMIT;
    //创建一个局部搜索器 built
    SlvaGenerator built = {.trellis=trellis, .limits=*limits, .next_sequence=1};
    SlvaStatus status;
    built.cells = allocate(&built, count*sizeof(*built.cells), &status);
    if (!built.cells) goto fail;

    status = reserve(&built, 1, length);// 为最小堆准备至少一个候选位置，同时确保预算中留有一条路径所需的 length 字节 bytes
    if (status != SLVA_OK) goto fail;

    uint8_t *bits = length ? allocate(&built, length, &status) : NULL;
    if (length && !bits) goto fail;

    //现在堆中只有一个候选节点，因为.free_length = length，代表整个合法路径集合。
    built.heap[0] = (SlvaCandidate){.bits=bits, .free_length=length};
    built.size = 1;
    BitVector path = {length,bits};

    //调用S12的viterbi_decode() 找到最小代价路径，并将其存储在 bits 中，完成初始化
    TrellisStatus decoded = viterbi_decode(trellis, built.cells, count, &path, &built.heap[0].cost);

    if (decoded == TRELLIS_NO_PATH) {
        release_path(&built, bits);
        built.size = 0;
    } else if (decoded != TRELLIS_OK) { status=from_trellis(decoded); goto fail; }
    *generator = built;

    //print_slva_generator(generator);
    return SLVA_OK;
fail:
    slva_free(&built);
    return status;
}

//从右向左逐个尝试“如果只把当前位置作为最后不同的位置，能否找到合法替代路径”
//返回当前节点能够划分出的非空子组数量，stop_at_one 只用于快速判断是否至少有一个子组。
static size_t alternative_count(const SlvaGenerator *g, const SlvaCandidate *node, int stop_at_one)
{
    const SyndromeTrellis *t = g->trellis;
    size_t state = t->target, count = 0;
    for (size_t j = t->length; j-- > 0;) {
        uint8_t bit = node->bits[j];
        size_t previous = state ^ (bit ? t->column_masks[j] : 0);
        size_t alternative = state ^ (bit ? 0 : t->column_masks[j]);// alternative prefix state
        if (j < node->free_length && isfinite(g->cells[j*t->state_count+alternative].metric)) {
            ++count;
            if (stop_at_one) return count;
        }
        state = previous;
    }
    return count;
}

//对一条已知完整路径重新计算代价并验证是否合法
static SlvaStatus score_path(const SyndromeTrellis *t, const uint8_t *bits, double *cost)
{
    size_t state = 0;
    double total = 0;
    for (size_t j = 0; j < t->length; ++j) {
        double edge;
        /*printf("第%zu层，bit=%u，state=%zu\n", j+1, bits[j], state);
        printf("计算下一状态和分支代价...\n");*/
        TrellisStatus status = syndrome_trellis_branch(t, j, state, bits[j], &state, &edge);
        if (status != TRELLIS_OK) return from_trellis(status);
        total += edge;
        if (!isfinite(total)) return SLVA_NUMERIC_ERROR;
    }
    if (state != t->target) return SLVA_INVALID_ARGUMENT;
    *cost = total;
    return SLVA_OK;
}

//把刚刚输出路径所在分组的剩余部分，按照最后不同位置划分成若干新子组，并为每个子组生成最轻代表。
static SlvaStatus expand_pending(SlvaGenerator *g)
{
    //printf("Expanding pending candidate...\n");
    print_slva_generator(g);
    if (!g->pending.bits) return SLVA_OK;
    const SyndromeTrellis *t = g->trellis;

    //检查可以产生几个子组，并预留资源
    size_t count = alternative_count(g, &g->pending, 0), live = g->size+1;
    if (live > g->limits.max_candidates || count > g->limits.max_candidates-live ||count > UINT64_MAX-g->next_sequence || (t->length && count > SIZE_MAX/t->length)) 
        return SLVA_RESOURCE_LIMIT;
    SlvaStatus status = reserve(g, g->size+count, count*t->length);
    if (status != SLVA_OK) 
        return status;
    
    //逐个寻找替代分支，从右到左沿原路径回溯。
    size_t state = t->target, built = 0, old_size = g->size;
    //printf("length = %zu\n",t->length);
    for (size_t j = t->length; j-- > 0;) {
        //printf("length = %zu\n",t->length);
        uint8_t bit = g->pending.bits[j];
        /*for(size_t i=0;i<t->length;i++){
            printf("%d ",g->pending.bits[i]);
        }
        printf("\n");*/
        size_t previous = state ^ (bit ? t->column_masks[j] : 0);
        size_t alternative = state ^ (bit ? 0 : t->column_masks[j]);
        printf("column_masks[%zu]=%zu\n",j,t->column_masks[j]);
        printf("state=%zu,previous=%zu,alternative=%zu\n",state,previous,alternative);
        
        printf("**free_length of the pending =%zu**\n",g->pending.free_length);
        
        if (j < g->pending.free_length && isfinite(g->cells[j*t->state_count+alternative].metric)) {
            uint8_t *bits = allocate(g, t->length, &status);
            if (!bits) goto fail;

            //翻转第j位，通过Viterbi前缀表补全新路径，生成最轻代表路径
            memcpy(bits, g->pending.bits, t->length);
            bits[j] ^= 1;
            size_t prefix_state = alternative;//把翻转后所需的前驱状态作为回溯起点
            //利用S12保存的前缀表，逐层恢复到达该状态的最轻前缀
            for (size_t layer = j; layer > 0; --layer) {
                const ViterbiCell *cell = &g->cells[layer*t->state_count+prefix_state];
                bits[layer-1] = cell->bit;
                prefix_state = cell->previous_state;
            }
            printf("Built new candidate with last differing position at %zu: ", j);
            printf("bits = ");
            for (size_t i = 0; i < t->length; ++i) printf("%d", bits[i]);
            printf("\n");
            printf("sequence = %" PRIu64 "\n", g->next_sequence+built-1);
            printf("free_length = %zu\n", j+1);
            

            //bits 就是新子组的最轻代表路径，计算其代价并放入堆中
            SlvaCandidate *child = &g->heap[old_size+built++];
            *child = (SlvaCandidate){.bits=bits, .free_length=j, .sequence=g->next_sequence+built-1};//生成一个候选分组，按.free_length = j定义（包含最小路径）
            status = score_path(t, bits, &child->cost);
            printf("cost = %f\n",child->cost);
            printf("----------------------------------\n");
            if (status != SLVA_OK) goto fail;
        }
        else {printf("No more candidates to build\n");}
        state = previous;
    }
    //逐个插入新生成得代表
    for (size_t i = 0; i < built; ++i) {
        SlvaCandidate child = g->heap[old_size+i]; //候选分组放入堆中
        push(g, child);
    }
    print_slva_heap(g->heap, g->size);

    g->next_sequence += built;
    release_path(g, g->pending.bits);
    g->pending = (SlvaCandidate){0};
    return SLVA_OK;
fail:
    for (size_t i = 0; i < built; ++i) release_path(g, g->heap[old_size+i].bits);
    return status;
}

/*确保下一条候选已经准备好
 *主要负责判断是否已经耗尽、是否达到路径输出上限，以及是否需要扩展 pending*/
static SlvaStatus ready(SlvaGenerator *g)
{
    if (!g || !g->trellis) 
        return SLVA_INVALID_ARGUMENT;
    if (!g->size && !g->pending.bits) 
        return SLVA_DONE;
    if (g->limits.path_limit && g->paths_returned >= g->limits.path_limit) 
        return SLVA_PATH_LIMIT;
    if (g->paths_returned == UINT64_MAX) 
        return SLVA_RESOURCE_LIMIT;
    SlvaStatus status = expand_pending(g);//必须先扩展，才能确定下一条真正的全局最优候选
    if (status != SLVA_OK) return status;
    return g->size ? SLVA_OK : SLVA_DONE;
}

//查看下一条候选的代价，但不输出路径，用于后续早停判断
SlvaStatus slva_peek(SlvaGenerator *g, double *cost)
{
    if (!cost) 
        return SLVA_INVALID_ARGUMENT;
    SlvaStatus status = ready(g);
    if (status == SLVA_OK) 
        *cost = g->heap[0].cost;
    return status;
}

//真正输出下一条候选路径
SlvaStatus slva_next(SlvaGenerator *g, BitVector *path, double *cost)
{
    if (!g || !g->trellis || !path || !cost || (g->trellis->length && !path->data)) 
        return SLVA_INVALID_ARGUMENT;
    if (path->length != g->trellis->length) 
        return SLVA_DIMENSION_MISMATCH;
    SlvaStatus status = ready(g);
    if (status != SLVA_OK) 
        return status;
    SlvaCandidate node = pop(g);
    if (path->length) 
        memcpy(path->data, node.bits, path->length);
    *cost = node.cost;
    ++g->paths_returned;
    if (alternative_count(g, &node, 1)) //输出最优分组代表，并保留剩余分组等待扩展
    else release_path(g, node.bits);

    /*printf("当前输出的候选路径为：\n");
    for (size_t i = 0; i < path->length; ++i)
        printf("%d ", path->data[i]);
    printf("\n");  
    printf("当前候选路径的代价为：%.2lf\n", *cost);*/
    

    return SLVA_OK;
}
