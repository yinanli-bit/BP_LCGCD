#include "bplcgcd/gf2.h"

#include <stdlib.h>
#include <string.h>

static unsigned checks, matrices;
#define CHECK(condition) do { \
    if (!(condition)) { \
        fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, #condition); \
        exit(EXIT_FAILURE); \
    } \
    ++checks; \
} while (0)

/* Independent integer arithmetic oracle, without elimination or library calls. */
static unsigned syndrome(const BinaryMatrix *h, unsigned mask)
{
    unsigned result = 0;
    for (size_t i = 0; i < h->rows; ++i) {
        unsigned sum = 0;
        for (size_t j = 0; j < h->cols; ++j) {
            sum += h->data[i * h->cols + j] * ((mask >> j) & 1u);
        }
        result |= (sum % 2) << i;
    }
    return result;
}

static unsigned reference_rank(const BinaryMatrix *h)
{
    unsigned seen[16] = {0}, count = 0, rank = 0;
    for (unsigned x = 0; x < (1u << h->cols); ++x) seen[syndrome(h, x)] = 1;
    for (unsigned i = 0; i < (1u << h->rows); ++i) count += seen[i];
    while ((1u << rank) < count) ++rank;
    CHECK((1u << rank) == count);
    return rank;
}

static void verify(const BinaryMatrix *h)
{
    BinarySystematic f = {0};
    BinaryMatrix uh = {0};
    BitVector e = {0}, es = {0}, restored = {0}, s = {0}, ts = {0}, rs = {0};
    uint8_t original[28];
    CHECK(h->rows <= 4 && h->cols <= 7);
    memcpy(original, h->data, h->rows * h->cols);
    CHECK(binary_matrix_systematize(h, &f) == GF2_OK);
    CHECK(f.rank == reference_rank(h));
    CHECK(reference_rank(&f.row_transform) == h->rows); /* U is invertible. */
    CHECK(binary_matrix_init(&uh, h->rows, h->cols) == GF2_OK);
    CHECK(binary_matrix_multiply(&f.row_transform, h, &uh) == GF2_OK);
    for (size_t j = 0; j < h->cols; ++j) {
        CHECK(f.sys_to_orig[j] < h->cols);
        CHECK(f.orig_to_sys[f.sys_to_orig[j]] == j);
        for (size_t i = 0; i < h->rows; ++i) {
            uint8_t value = f.reduced.data[i * h->cols + j];
            CHECK(value == uh.data[i * h->cols + f.sys_to_orig[j]]);
            if (j < f.rank) CHECK(value == (i == j));
            if (i >= f.rank) CHECK(value == 0);
        }
    }
    CHECK(bit_vector_init(&e, h->cols) == GF2_OK);
    CHECK(bit_vector_init(&es, h->cols) == GF2_OK);
    CHECK(bit_vector_init(&restored, h->cols) == GF2_OK);
    CHECK(bit_vector_init(&s, h->rows) == GF2_OK);
    CHECK(bit_vector_init(&ts, h->rows) == GF2_OK);
    CHECK(bit_vector_init(&rs, h->rows) == GF2_OK);
    for (unsigned x = 0; x < (1u << h->cols); ++x) {
        for (size_t j = 0; j < h->cols; ++j) e.data[j] = (uint8_t)((x >> j) & 1u);
        CHECK(binary_systematic_to_system(&f, &e, &es) == GF2_OK);
        CHECK(binary_systematic_to_original(&f, &es, &restored) == GF2_OK);
        CHECK(memcmp(e.data, restored.data, h->cols) == 0);
        for (size_t j = 0; j < h->cols; ++j) CHECK(es.data[j] == e.data[f.sys_to_orig[j]]);
        CHECK(binary_matrix_multiply_vector(h, &e, &s) == GF2_OK);
        CHECK(binary_systematic_transform_syndrome(&f, &s, &ts) == GF2_OK);
        CHECK(binary_matrix_multiply_vector(&f.reduced, &es, &rs) == GF2_OK);
        CHECK(memcmp(ts.data, rs.data, h->rows) == 0);
        for (size_t i = 0; i < h->rows; ++i) CHECK(s.data[i] == ((syndrome(h, x) >> i) & 1u));
    }
    for (unsigned target = 0; target < (1u << h->rows); ++target) {
        int exists = 0;
        for (unsigned x = 0; x < (1u << h->cols); ++x) {
            if (syndrome(h, x) == target) exists = 1;
        }
        for (size_t i = 0; i < h->rows; ++i) s.data[i] = (uint8_t)((target >> i) & 1u);
        CHECK(binary_systematic_transform_syndrome(&f, &s, &ts) ==
              (exists ? GF2_OK : GF2_INCONSISTENT_SYNDROME));
        for (size_t i = 0; i < h->rows; ++i) {
            CHECK(ts.data[i] == ((syndrome(&f.row_transform, target) >> i) & 1u));
        }
    }
    CHECK(memcmp(h->data, original, h->rows * h->cols) == 0);
    binary_matrix_free(&uh);
    bit_vector_free(&e); bit_vector_free(&es); bit_vector_free(&restored);
    bit_vector_free(&s); bit_vector_free(&ts); bit_vector_free(&rs);
    binary_systematic_free(&f);
    binary_systematic_free(&f);
    CHECK(!f.rank && !f.reduced.data && !f.row_transform.data &&
          !f.sys_to_orig && !f.orig_to_sys);
    ++matrices;
}

static void test_errors(void)
{
    BinaryMatrix h = {0}, empty = {0};
    BinarySystematic f = {0}, absent = {0};
    BitVector e = {0}, out = {0}, s = {0}, ts = {0};
    CHECK(binary_matrix_init(&h, 2, 3) == GF2_OK);
    CHECK(binary_matrix_systematize(NULL, &f) == GF2_INVALID_ARGUMENT);
    CHECK(binary_matrix_systematize(&h, NULL) == GF2_INVALID_ARGUMENT);
    CHECK(binary_matrix_systematize(&empty, &f) == GF2_INVALID_ARGUMENT);
    h.data[5] = 2;
    CHECK(binary_matrix_systematize(&h, &f) == GF2_INVALID_ARGUMENT);
    CHECK(!f.reduced.data && !f.row_transform.data && !f.sys_to_orig);
    h.data[5] = 0;
    CHECK(binary_matrix_systematize(&h, &f) == GF2_OK);
    uint8_t *saved = f.reduced.data;
    CHECK(binary_matrix_systematize(&h, &f) == GF2_INVALID_ARGUMENT);
    CHECK(saved == f.reduced.data);
    CHECK(bit_vector_init(&e, 3) == GF2_OK);
    CHECK(bit_vector_init(&out, 3) == GF2_OK);
    CHECK(bit_vector_init(&s, 2) == GF2_OK);
    CHECK(bit_vector_init(&ts, 2) == GF2_OK);
    memset(out.data, 1, 3);
    e.data[2] = 2;
    CHECK(binary_systematic_to_system(&f, &e, &out) == GF2_INVALID_ARGUMENT);
    CHECK(binary_systematic_to_original(&f, &e, &out) == GF2_INVALID_ARGUMENT);
    CHECK(out.data[0] == 1 && out.data[1] == 1 && out.data[2] == 1);
    e.data[2] = 0;
    CHECK(binary_systematic_to_system(&f, &e, &e) == GF2_INVALID_ARGUMENT);
    CHECK(binary_systematic_to_system(&absent, &e, &out) == GF2_INVALID_ARGUMENT);
    CHECK(binary_systematic_to_original(&f, NULL, &out) == GF2_INVALID_ARGUMENT);
    CHECK(binary_systematic_to_original(&f, &e, NULL) == GF2_INVALID_ARGUMENT);
    CHECK(binary_systematic_to_system(&f, &s, &out) == GF2_DIMENSION_MISMATCH);
    CHECK(binary_systematic_to_original(&f, &e, &ts) == GF2_DIMENSION_MISMATCH);
    CHECK(binary_systematic_transform_syndrome(&absent, &s, &ts) == GF2_INVALID_ARGUMENT);
    CHECK(binary_systematic_transform_syndrome(&f, &s, &s) == GF2_INVALID_ARGUMENT);
    s.data[1] = 2;
    CHECK(binary_systematic_transform_syndrome(&f, &s, &ts) == GF2_INVALID_ARGUMENT);
    CHECK(ts.data[0] == 0 && ts.data[1] == 0);
    s.data[1] = 1;
    CHECK(binary_systematic_transform_syndrome(&f, &e, &ts) == GF2_DIMENSION_MISMATCH);
    CHECK(binary_systematic_transform_syndrome(&f, &s, &ts) == GF2_INCONSISTENT_SYNDROME);
    CHECK(ts.data[0] == 0 && ts.data[1] == 1); /* rank=0 retains diagnostics. */
    bit_vector_free(&e); bit_vector_free(&out); bit_vector_free(&s); bit_vector_free(&ts);
    binary_matrix_free(&h);
    binary_systematic_free(&f);
    binary_systematic_free(NULL);
}

int main(void)
{
    test_errors();
    for (size_t rows = 1; rows <= 3; ++rows) {
        for (size_t cols = 1; cols <= 3; ++cols) {
            BinaryMatrix h = {0};
            CHECK(binary_matrix_init(&h, rows, cols) == GF2_OK);
            for (unsigned mask = 0; mask < (1u << (rows * cols)); ++mask) {
                for (size_t i = 0; i < rows * cols; ++i) h.data[i] = (uint8_t)((mask >> i) & 1u);
                verify(&h);
            }
            binary_matrix_free(&h);
        }
    }
    BinaryMatrix h = {0};
    CHECK(binary_matrix_init(&h, 3, 4) == GF2_OK);
    const uint8_t base[] = {0,1,1,0, 0,0,1,1, 0,1,0,1};
    unsigned permutations = 0;
    for (unsigned code = 0; code < 256; ++code) {
        unsigned seen = 0;
        for (size_t j = 0; j < 4; ++j) seen |= 1u << ((code >> (2 * j)) & 3u);
        if (seen != 15) continue; /* Exactly the 24 permutations of four columns. */
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 4; ++j) h.data[i * 4 + j] = base[i * 4 + ((code >> (2 * j)) & 3u)];
        }
        verify(&h);
        ++permutations;
    }
    CHECK(permutations == 24);
    memcpy(h.data, base, sizeof(base));
    BinarySystematic example = {0};
    CHECK(binary_matrix_systematize(&h, &example) == GF2_OK);
    const uint8_t reduced[] = {1,0,0,1, 0,1,0,1, 0,0,0,0};
    CHECK(example.rank == 2 && memcmp(example.reduced.data, reduced, sizeof(reduced)) == 0);
    CHECK(example.sys_to_orig[0] == 1 && example.sys_to_orig[1] == 2 &&
          example.sys_to_orig[2] == 0 && example.sys_to_orig[3] == 3);
    puts("Hand example: rank=2, K=2, sys_to_orig=(1,2,0,3); redundant row retained");
    binary_systematic_free(&example);
    binary_matrix_free(&h);
    CHECK(binary_matrix_init(&h, 3, 7) == GF2_OK);
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 7; ++j) h.data[i * 7 + j] = (uint8_t)(((j + 1) >> i) & 1u);
    }
    verify(&h); /* Hamming [7,4], all 128 vectors and all eight syndromes. */
    binary_matrix_free(&h);
    CHECK(binary_matrix_init(&h, 4, 2) == GF2_OK);
    const uint8_t tall[] = {0,1, 1,0, 1,1, 0,0};
    memcpy(h.data, tall, sizeof(tall));
    verify(&h);
    binary_matrix_free(&h);
    printf("Systematic checks passed (%u checks, %u matrices, 24 column permutations)\n", checks, matrices);
    return EXIT_SUCCESS;
}
