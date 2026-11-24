#include "bplcgcd/slva.h"
#include "bplcgcd/tep.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static unsigned checks, graphs;
#define CHECK(condition) do { \
    if (!(condition)) { \
        fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, #condition); \
        exit(EXIT_FAILURE); \
    } \
    ++checks; \
} while (0)

typedef struct { unsigned bits; double cost; } OraclePath;
static const SlvaLimits generous = {4096, 1024*1024, 0};

static int compare_cost(const void *a, const void *b)
{
    const OraclePath *x=a, *y=b;
    if (x->cost != y->cost) return x->cost < y->cost ? -1 : 1;
    return x->bits < y->bits ? -1 : x->bits != y->bits;
}

static unsigned pack(const BitVector *path)
{
    unsigned bits = 0;
    for (size_t j = 0; j < path->length; ++j) {
        CHECK(path->data[j] <= 1);
        bits |= (unsigned)path->data[j] << j;
    }
    return bits;
}

/* Independent Eq.9 vector arithmetic: no trellis masks/branch or SLVA code. */
static int local_valid(const LocalConstraints *local, unsigned bits)
{
    size_t delta=local->p_a.rows, k=local->p_a.cols;
    for (size_t i = 0; i < delta; ++i) {
        unsigned parity = (bits >> i)&1u;
        for (size_t j = 0; j < k; ++j) parity ^= local->p_a.data[i*k+j] & ((bits >> (delta+j))&1u);
        if (parity != local->s_a.data[i]) return 0;
    }
    return 1;
}

static void verify(const BinarySystematic *form, const double *weights, const BitVector *residual, size_t delta)
{
    size_t m=form->rank, k=form->reduced.cols-m, length=delta+k;
    size_t order[8], masks[8], columns[8];
    uint8_t pa[32], pac[32], sa[4], sac[4], scratch_data[4], bits[8], replay_bits[8];
    LocalConstraints local = {.basis_order=m ? order : NULL,
        .p_a={delta,k,delta && k ? pa : NULL}, .p_ac={m-delta,k,m>delta && k ? pac : NULL},
        .s_a={delta,delta ? sa : NULL}, .s_ac={m-delta,m>delta ? sac : NULL}};
    BitVector scratch = {form->reduced.rows,scratch_data};
    CHECK(local_constraints_prepare(form,weights,form->reduced.cols,residual,delta,&local,&scratch) == GF2_OK);
    SyndromeTrellis t = {.length=length,.column_masks=masks,.original_columns=columns};
    CHECK(syndrome_trellis_prepare(form,&local,weights,form->reduced.cols,16,&t) == TRELLIS_OK);
    OraclePath oracle[256];
    unsigned count=0, emitted[256];
    uint8_t seen[256] = {0};
    for (unsigned word = 0; word < (1u << length); ++word) {
        if (!local_valid(&local,word)) continue;
        double cost = 0;
        for (size_t j = 0; j < length; ++j) {
            size_t original = j < delta ? order[j] : form->sys_to_orig[m+j-delta];
            cost += ((word >> j)&1u)*weights[original];
        }
        oracle[count++] = (OraclePath){word,cost};
    }
    CHECK(count == (1u << k));
    qsort(oracle,count,sizeof(*oracle),compare_cost);
    SlvaGenerator g = {0};
    CHECK(slva_init(&g,&t,&generous) == SLVA_OK);
    BitVector path = {length,length ? bits : NULL};
    ViterbiCell cells[144];
    BitVector first = {length,length ? replay_bits : NULL};
    double first_cost;
    CHECK(viterbi_decode(&t,cells,144,&first,&first_cost) == TRELLIS_OK);
    for (unsigned rank = 0; rank < count; ++rank) {
        double cost=-1, peek=-2, again=-3;
        CHECK(slva_peek(&g,&peek) == SLVA_OK && g.paths_returned == rank);
        CHECK(slva_peek(&g,&again) == SLVA_OK && again == peek && g.paths_returned == rank);
        CHECK(slva_next(&g,&path,&cost) == SLVA_OK && cost == peek);
        CHECK(cost == oracle[rank].cost); /* Equal-cost paths may use any stable order. */
        unsigned word = pack(&path);
        CHECK(!seen[word] && local_valid(&local,word));
        seen[word]=1; emitted[rank]=word;
        double direct = 0;
        for (size_t j = 0; j < length; ++j) direct += bits[j]*weights[columns[j]];
        CHECK(cost == direct && g.paths_returned == rank+1);
        CHECK(g.bytes_used <= generous.max_bytes && g.size+(g.pending.bits != NULL) <= generous.max_candidates);
        if (!rank) CHECK(cost == first_cost && !memcmp(bits,replay_bits,length));
    }
    for (unsigned i = 0; i < count; ++i) CHECK(seen[oracle[i].bits]);
    memset(bits,9,sizeof(bits));
    double cost=99;
    CHECK(slva_next(&g,&path,&cost) == SLVA_DONE && cost == 99 && bits[0] == 9);
    CHECK(slva_peek(&g,&cost) == SLVA_DONE && cost == 99);
    CHECK(slva_next(&g,&path,&cost) == SLVA_DONE);
    slva_free(&g);
    CHECK(!g.trellis && !g.cells && !g.heap && !g.bytes_used);
    /* Reinitialization and prefix limits must reproduce the very same list. */
    const uint64_t limits[] = {0,1,3,17,count};
    for (size_t trial = 0; trial < sizeof(limits)/sizeof(*limits); ++trial) {
        SlvaLimits bounded = generous;
        bounded.path_limit=limits[trial];
        CHECK(slva_init(&g,&t,&bounded) == SLVA_OK);
        unsigned take = bounded.path_limit && bounded.path_limit<count ? (unsigned)bounded.path_limit : count;
        for (unsigned rank = 0; rank < take; ++rank) {
            CHECK(slva_next(&g,&path,&cost) == SLVA_OK);
            CHECK(pack(&path) == emitted[rank] && cost == oracle[rank].cost);
        }
        SlvaStatus stop = take<count ? SLVA_PATH_LIMIT : SLVA_DONE;
        CHECK(slva_next(&g,&path,&cost) == stop && g.paths_returned == take);
        CHECK(slva_peek(&g,&cost) == stop);
        slva_free(&g);
    }
    ++graphs;
}

static void test_manual(void)
{
    uint8_t data[] = {0,1,1,0,1,0, 0,0,1,1,0,1, 0,1,0,1,1,1};
    BinaryMatrix h = {3,6,data};
    BinarySystematic form = {0};
    CHECK(binary_matrix_systematize(&h,&form) == GF2_OK);
    double weights[] = {1000,3,6,9,8,7};
    uint8_t residual_data[] = {1,1,0};
    BitVector residual = {3,residual_data};
    for (size_t delta = 0; delta <= form.rank; ++delta) verify(&form,weights,&residual,delta);
    size_t masks[] = {1,2,0,3,2,3}, columns[] = {2,1,0,3,4,5};
    SyndromeTrellis t = {.length=6,.delta=2,.state_count=4,.target=1,
        .column_masks=masks,.original_columns=columns,.flip_weights=weights};
    SlvaGenerator g = {0};
    CHECK(slva_init(&g,&t,&generous) == SLVA_OK);
    uint8_t bits[6];
    BitVector path = {6,bits};
    const uint8_t expected[3][6] = {{1,0,0,0,0,0},{0,1,0,0,0,1},{0,1,0,1,0,0}};
    const double costs[] = {6,10,12};
    for (size_t rank = 0; rank < 3; ++rank) {
        double cost;
        CHECK(slva_next(&g,&path,&cost) == SLVA_OK);
        CHECK(cost == costs[rank] && !memcmp(bits,expected[rank],6));
        printf("SLVA hand rank=%zu, bits=",rank+1);
        for (size_t j = 0; j < 6; ++j) printf("%u",bits[j]);
        printf(", local cost=%.0f\n",cost);
    }
    slva_free(&g);
    binary_systematic_free(&form);
}

static void test_limits_and_empty(void)
{
    size_t masks[] = {0,0}, columns[] = {0,1};
    double weights[] = {1,2}, cost=99;
    SyndromeTrellis t = {.length=2,.delta=0,.state_count=1,.target=0,
        .column_masks=masks,.original_columns=columns,.flip_weights=weights};
    uint8_t bits[] = {9,9};
    BitVector path = {2,bits};
    SlvaGenerator g = {0};
    CHECK(slva_init(NULL,&t,&generous) == SLVA_INVALID_ARGUMENT);
    CHECK(slva_init(&g,NULL,&generous) == SLVA_INVALID_ARGUMENT);
    CHECK(slva_init(&g,&t,NULL) == SLVA_INVALID_ARGUMENT);
    SlvaLimits budget = generous;
    budget.max_candidates=0;
    CHECK(slva_init(&g,&t,&budget) == SLVA_INVALID_ARGUMENT);
    budget=generous; budget.max_bytes=1;
    CHECK(slva_init(&g,&t,&budget) == SLVA_RESOURCE_LIMIT && !g.trellis);
    for (unsigned trial = 0; trial < 3; ++trial) {
        budget=generous;
        if (trial<2) budget.max_candidates=trial+1; /* Expansion also counts the pending parent. */
        else budget.max_bytes=3*sizeof(ViterbiCell)+sizeof(SlvaCandidate)+2;
        CHECK(slva_init(&g,&t,&budget) == SLVA_OK);
        CHECK(slva_next(&g,&path,&cost) == SLVA_OK && cost == 0);
        CHECK(g.paths_returned == 1 && g.pending.bits);
        memset(bits,9,2); cost=99;
        CHECK(slva_peek(&g,&cost) == SLVA_RESOURCE_LIMIT && cost == 99);
        CHECK(slva_next(&g,&path,&cost) == SLVA_RESOURCE_LIMIT && cost == 99 && bits[0] == 9);
        CHECK(g.paths_returned == 1 && g.pending.bits && g.size == 0);
        CHECK(g.bytes_used <= budget.max_bytes);
        slva_free(&g);
    }
    budget=generous; budget.max_candidates=1; budget.path_limit=1;
    CHECK(slva_init(&g,&t,&budget) == SLVA_OK);
    CHECK(slva_next(&g,&path,&cost) == SLVA_OK && cost == 0);
    CHECK(slva_next(&g,&path,&cost) == SLVA_PATH_LIMIT && g.paths_returned == 1);
    slva_free(&g); /* First-path use never pays for a second-path expansion. */
    CHECK(slva_init(&g,&t,&generous) == SLVA_OK);
    CHECK(slva_next(NULL,&path,&cost) == SLVA_INVALID_ARGUMENT);
    CHECK(slva_next(&g,NULL,&cost) == SLVA_INVALID_ARGUMENT);
    CHECK(slva_next(&g,&path,NULL) == SLVA_INVALID_ARGUMENT);
    CHECK(slva_peek(&g,NULL) == SLVA_INVALID_ARGUMENT);
    path.length=1;
    CHECK(slva_next(&g,&path,&cost) == SLVA_DIMENSION_MISMATCH && g.paths_returned == 0);
    path.length=2;
    CHECK(slva_next(&g,&path,&cost) == SLVA_OK && cost == 0);
    slva_free(&g);
    t.delta=1; t.state_count=2; t.target=1; /* No path: two zero columns. */
    CHECK(slva_init(&g,&t,&generous) == SLVA_OK);
    cost=99;
    CHECK(slva_next(&g,&path,&cost) == SLVA_DONE && cost == 99 && g.paths_returned == 0);
    slva_free(&g);
    t.length=path.length=0; t.state_count=1; t.target=0; path.data=NULL;
    budget=generous; budget.path_limit=1;
    CHECK(slva_init(&g,&t,&budget) == SLVA_OK);
    CHECK(slva_next(&g,&path,&cost) == SLVA_OK && cost == 0);
    CHECK(slva_next(&g,&path,&cost) == SLVA_DONE && g.paths_returned == 1);
    slva_free(&g);
    t.length=SIZE_MAX;
    CHECK(slva_init(&g,&t,&generous) == SLVA_RESOURCE_LIMIT && !g.cells);
    t.length=1; t.state_count=(size_t)1 << (sizeof(size_t)*CHAR_BIT-1);
    CHECK(slva_init(&g,&t,&generous) == SLVA_RESOURCE_LIMIT && !g.cells);
    t.length=0;
    CHECK(slva_init(&g,&t,&generous) == SLVA_RESOURCE_LIMIT && !g.cells);
    t.length=2; t.state_count=2; t.target=0; masks[0]=masks[1]=1;
    weights[0]=weights[1]=DBL_MAX;
    CHECK(slva_init(&g,&t,&generous) == SLVA_NUMERIC_ERROR && !g.cells && !g.heap);
    slva_free(&g); slva_free(&g); slva_free(NULL);
}

static void test_long_prefix(void)
{
    size_t masks[70] = {0}, columns[70];
    double weights[70];
    for (size_t j = 0; j < 70; ++j) { columns[j]=j; weights[j]=(double)(70-j); }
    SyndromeTrellis t = {.length=70,.state_count=1,.column_masks=masks,
        .original_columns=columns,.flip_weights=weights};
    SlvaLimits limits = generous;
    limits.path_limit=25;
    SlvaGenerator g = {0};
    TepGenerator tep = {0};
    CHECK(slva_init(&g,&t,&limits) == SLVA_OK);
    CHECK(tep_generator_init(&tep,weights,70) == TEP_OK);
    uint8_t bits[70], tep_bits[70], history[25][70];
    BitVector path = {70,bits}, tep_path = {70,tep_bits};
    for (size_t rank = 0; rank < 25; ++rank) {
        double cost, tep_cost, direct=0;
        CHECK(slva_next(&g,&path,&cost) == SLVA_OK);
        CHECK(tep_generator_next(&tep,&tep_path,&tep_cost) == TEP_OK && cost == tep_cost);
        for (size_t j = 0; j < 70; ++j) { CHECK(bits[j] <= 1); direct += bits[j]*weights[j]; }
        CHECK(cost == direct && g.bytes_used <= limits.max_bytes);
        for (size_t previous = 0; previous < rank; ++previous) CHECK(memcmp(bits,history[previous],70));
        memcpy(history[rank],bits,70);
    }
    double cost=99;
    CHECK(slva_next(&g,&path,&cost) == SLVA_PATH_LIMIT && cost == 99);
    slva_free(&g); tep_generator_free(&tep);
    CHECK(slva_init(&g,&t,&limits) == SLVA_OK);
    for (size_t rank = 0; rank < 25; ++rank) {
        CHECK(slva_next(&g,&path,&cost) == SLVA_OK && !memcmp(bits,history[rank],70));
    }
    slva_free(&g);
    puts("70-bit delta=0: first 25 SLVA costs match TEP; deterministic, unique, bounded prefix");
}

int main(void)
{
    test_manual();
    test_limits_and_empty();
    test_long_prefix();
    uint32_t seed=20261123u;
    for (unsigned trial = 0; trial < 48; ++trial) {
        uint8_t data[32] = {0}, residual_data[4] = {0};
        BinaryMatrix h = {1+trial%4,1+trial%8,data};
        if (trial == 1) { h.rows=h.cols=3; data[0]=data[4]=data[8]=1; } /* K=0. */
        else if (trial%8) {
            for (size_t j = 0; j < h.rows*h.cols; ++j) { seed=1664525u*seed+1013904223u; data[j]=(uint8_t)(seed>>31); }
        } /* Zero H includes rank=0 and all labeled paths. */
        for (size_t i = 0; i < h.rows; ++i) {
            for (size_t j = 0; j < h.cols; ++j) residual_data[i] ^= data[i*h.cols+j] & (uint8_t)((37u>>j)&1u);
        }
        BinarySystematic form = {0};
        CHECK(binary_matrix_systematize(&h,&form) == GF2_OK);
        double weights[8];
        for (size_t j = 0; j < h.cols; ++j) {
            const double decimal[] = {0,1e15,0.1,0.2,0.3,1e-14,0.4,1e-8};
            if (trial%4 == 0) weights[j]=0;
            else if (trial%4 == 1) weights[j]=1;
            else if (trial%4 == 2) weights[j]=(double)((3*j+trial)%7)/4;
            else weights[j]=decimal[(j+trial)%8];
        }
        BitVector residual = {h.rows,residual_data};
        for (size_t delta = 0; delta <= form.rank; ++delta) verify(&form,weights,&residual,delta);
        binary_systematic_free(&form);
    }
    printf("SLVA checks passed (%u checks, %u exhaustive graphs; complete lists, ties, prefixes and budgets)\n",checks,graphs);
    return EXIT_SUCCESS;
}
