#include "bplcgcd/gf2.h"

#include <stdlib.h>
#include <string.h>

static unsigned checks;
#define CHECK(condition) do { \
    if (!(condition)) { \
        fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, #condition); \
        exit(EXIT_FAILURE); \
    } \
    ++checks; \
} while (0)

static Gf2Status load_text(const char *text, BinaryMatrix *matrix)
{
    FILE *stream = tmpfile();
    CHECK(stream != NULL);
    CHECK(fputs(text, stream) >= 0);
    CHECK(fflush(stream) == 0);
    rewind(stream);
    Gf2Status status = binary_matrix_load(stream, matrix);
    CHECK(fclose(stream) == 0);
    return status;
}

static void test_storage(void)
{
    BinaryMatrix h = {0};
    BitVector v = {0};
    CHECK(binary_matrix_init(NULL, 1, 1) == GF2_INVALID_ARGUMENT);
    CHECK(binary_matrix_init(&h, 0, 1) == GF2_INVALID_ARGUMENT);
    CHECK(binary_matrix_init(&h, 1, 0) == GF2_INVALID_ARGUMENT);
    CHECK(binary_matrix_init(&h, SIZE_MAX, 2) == GF2_INVALID_ARGUMENT);
    CHECK(h.data == NULL && h.rows == 0 && h.cols == 0);
    CHECK(bit_vector_init(NULL, 1) == GF2_INVALID_ARGUMENT);
    CHECK(bit_vector_init(&v, 0) == GF2_INVALID_ARGUMENT);
    CHECK(binary_matrix_init(&h, 2, 9) == GF2_OK);
    CHECK(bit_vector_init(&v, 9) == GF2_OK);
    for (size_t i = 0; i < 18; ++i) CHECK(h.data[i] == 0);
    for (size_t i = 0; i < 9; ++i) CHECK(v.data[i] == 0);
    uint8_t *original = h.data;
    CHECK(binary_matrix_init(&h, 1, 1) == GF2_INVALID_ARGUMENT);
    CHECK(h.data == original && h.rows == 2 && h.cols == 9);
    CHECK(bit_vector_init(&v, 1) == GF2_INVALID_ARGUMENT);
    binary_matrix_free(&h);
    binary_matrix_free(&h);
    binary_matrix_free(NULL);
    bit_vector_free(&v);
    bit_vector_free(&v);
    bit_vector_free(NULL);
    CHECK(h.data == NULL && h.rows == 0 && h.cols == 0);
    CHECK(v.data == NULL && v.length == 0);
}

static void test_load(void)
{
    const char *invalid[] = {
        "", " \t\r\n", "x 1 0", "0 1", "1 0", "-1 1 0", "+1 1 0",
        "1.0 1 0", "1 -1 0", "1 1", "2 2 0 1 1", "1 1 2", "1 1 -1",
        "1 1 x", "1 1 01", "1 1 1x", "1 1 0 1", "1 1 0 # comment"
    };
    BinaryMatrix h = {0};
    CHECK(binary_matrix_load(NULL, &h) == GF2_INVALID_ARGUMENT);
    CHECK(load_text("1 1 0", NULL) == GF2_INVALID_ARGUMENT);
    for (size_t i = 0; i < sizeof(invalid) / sizeof(invalid[0]); ++i) {
        CHECK(load_text(invalid[i], &h) == GF2_FORMAT_ERROR);
        CHECK(h.data == NULL && h.rows == 0 && h.cols == 0);
    }
    char overflow[128];
    CHECK(snprintf(overflow, sizeof(overflow), "%zu0 1 0", SIZE_MAX) > 0);
    CHECK(load_text(overflow, &h) == GF2_FORMAT_ERROR);
    CHECK(snprintf(overflow, sizeof(overflow), "%zu 2", SIZE_MAX) > 0);
    CHECK(load_text(overflow, &h) == GF2_FORMAT_ERROR);
    CHECK(h.data == NULL && h.rows == 0 && h.cols == 0);
    CHECK(load_text(" \t2 3\r\n1 0 1\r\n0\t1 0 \n", &h) == GF2_OK);
    const uint8_t expected[] = {1, 0, 1, 0, 1, 0};
    CHECK(h.rows == 2 && h.cols == 3);
    CHECK(memcmp(h.data, expected, sizeof(expected)) == 0);
    uint8_t *original = h.data;
    CHECK(load_text("1 1 0", &h) == GF2_INVALID_ARGUMENT);
    CHECK(h.data == original && memcmp(h.data, expected, sizeof(expected)) == 0);
    binary_matrix_free(&h);
    CHECK(load_text("1 1 1", &h) == GF2_OK); /* EOF without a final newline. */
    CHECK(h.data[0] == 1);
    binary_matrix_free(&h);
}

static void test_syndrome(void)
{
    BinaryMatrix h = {0};
    BitVector e = {0}, s = {0};
    CHECK(load_text("3 4\n1 0 1 1\n0 1 1 0\n1 1 0 1\n", &h) == GF2_OK);
    CHECK(bit_vector_init(&e, 4) == GF2_OK);
    CHECK(bit_vector_init(&s, 3) == GF2_OK);
    e.data[0] = e.data[2] = 1;
    CHECK(binary_matrix_multiply_vector(&h, &e, &s) == GF2_OK);
    CHECK(s.data[0] == 0 && s.data[1] == 1 && s.data[2] == 1);
    puts("Hand example: e=(1,0,1,0), H e^T=(0,1,1)");
    /* Independent integer bitmask oracle: bit j denotes column j. */
    const unsigned row_masks[] = {13, 6, 11};
    for (unsigned mask = 0; mask < 16; ++mask) {
        for (size_t j = 0; j < 4; ++j) e.data[j] = (uint8_t)((mask >> j) & 1u);
        CHECK(binary_matrix_multiply_vector(&h, &e, &s) == GF2_OK);
        for (size_t i = 0; i < 3; ++i) {
            unsigned value = mask & row_masks[i], sum = 0;
            while (value) { sum += value % 2; value /= 2; }
            CHECK(s.data[i] == sum % 2);
        }
    }
    /* Invalid inputs must leave even the first output entry unchanged. */
    const uint8_t sentinel[] = {1, 1, 1};
    memcpy(s.data, sentinel, sizeof(sentinel));
    h.data[11] = 2;
    CHECK(binary_matrix_multiply_vector(&h, &e, &s) == GF2_INVALID_ARGUMENT);
    CHECK(memcmp(s.data, sentinel, sizeof(sentinel)) == 0);
    h.data[11] = 1;
    e.data[3] = 2;
    CHECK(binary_matrix_multiply_vector(&h, &e, &s) == GF2_INVALID_ARGUMENT);
    CHECK(memcmp(s.data, sentinel, sizeof(sentinel)) == 0);
    e.data[3] = 1;
    CHECK(binary_matrix_multiply_vector(NULL, &e, &s) == GF2_INVALID_ARGUMENT);
    CHECK(binary_matrix_multiply_vector(&h, NULL, &s) == GF2_INVALID_ARGUMENT);
    CHECK(binary_matrix_multiply_vector(&h, &e, NULL) == GF2_INVALID_ARGUMENT);
    BitVector wrong = {0};
    CHECK(bit_vector_init(&wrong, 2) == GF2_OK);
    CHECK(binary_matrix_multiply_vector(&h, &wrong, &s) == GF2_DIMENSION_MISMATCH);
    CHECK(binary_matrix_multiply_vector(&h, &e, &wrong) == GF2_DIMENSION_MISMATCH);
    CHECK(memcmp(s.data, sentinel, sizeof(sentinel)) == 0);
    bit_vector_free(&wrong);
    binary_matrix_free(&h);
    bit_vector_free(&e);
    CHECK(binary_matrix_multiply_vector(&h, &e, &s) == GF2_INVALID_ARGUMENT);
    bit_vector_free(&s);
}

static void test_boundaries(void)
{
    BinaryMatrix h = {0};
    BitVector e = {0}, s = {0};
    CHECK(binary_matrix_init(&h, 1, 1) == GF2_OK);
    CHECK(bit_vector_init(&e, 1) == GF2_OK);
    CHECK(bit_vector_init(&s, 1) == GF2_OK);
    for (uint8_t a = 0; a < 2; ++a) {
        for (uint8_t b = 0; b < 2; ++b) {
            h.data[0] = a; e.data[0] = b;
            CHECK(binary_matrix_multiply_vector(&h, &e, &s) == GF2_OK);
            CHECK(s.data[0] == a * b);
        }
    }
    CHECK(binary_matrix_multiply_vector(&h, &e, &e) == GF2_INVALID_ARGUMENT);
    BitVector alias = {1, h.data};
    CHECK(binary_matrix_multiply_vector(&h, &e, &alias) == GF2_INVALID_ARGUMENT);
    binary_matrix_free(&h);
    bit_vector_free(&e);
    CHECK(binary_matrix_init(&h, 2, 9) == GF2_OK);
    CHECK(bit_vector_init(&e, 9) == GF2_OK);
    bit_vector_free(&s);
    CHECK(bit_vector_init(&s, 2) == GF2_OK);
    h.data[0] = h.data[8] = h.data[9 + 8] = 1;
    e.data[0] = e.data[8] = 1;
    CHECK(binary_matrix_multiply_vector(&h, &e, &s) == GF2_OK);
    CHECK(s.data[0] == 0 && s.data[1] == 1); /* Column 8 is a separate byte. */
    binary_matrix_free(&h);
    bit_vector_free(&e);
    bit_vector_free(&s);
}

static void test_matrix_product(void)
{
    BinaryMatrix a = {0}, b = {0}, c = {0}, empty = {0};
    CHECK(load_text("2 3 1 0 1 0 1 1", &a) == GF2_OK);
    CHECK(load_text("3 2 1 1 1 0 0 1", &b) == GF2_OK);
    CHECK(binary_matrix_init(&c, 2, 2) == GF2_OK);
    CHECK(binary_matrix_multiply(&a, &b, &c) == GF2_OK);
    const uint8_t expected[] = {1, 0, 1, 1};
    CHECK(memcmp(c.data, expected, sizeof(expected)) == 0);
    b.data[5] = 2;
    CHECK(binary_matrix_multiply(&a, &b, &c) == GF2_INVALID_ARGUMENT);
    b.data[5] = 1;
    a.data[5] = 2;
    CHECK(binary_matrix_multiply(&a, &b, &c) == GF2_INVALID_ARGUMENT);
    a.data[5] = 1;
    b.rows = 2;
    CHECK(binary_matrix_multiply(&a, &b, &c) == GF2_DIMENSION_MISMATCH);
    b.rows = 3;
    c.rows = 1;
    CHECK(binary_matrix_multiply(&a, &b, &c) == GF2_DIMENSION_MISMATCH);
    c.rows = 2;
    CHECK(binary_matrix_multiply(NULL, &b, &c) == GF2_INVALID_ARGUMENT);
    CHECK(binary_matrix_multiply(&a, NULL, &c) == GF2_INVALID_ARGUMENT);
    CHECK(binary_matrix_multiply(&a, &b, NULL) == GF2_INVALID_ARGUMENT);
    CHECK(binary_matrix_multiply(&empty, &b, &c) == GF2_INVALID_ARGUMENT);
    CHECK(binary_matrix_multiply(&a, &b, &empty) == GF2_INVALID_ARGUMENT);
    CHECK(binary_matrix_multiply(&a, &b, &a) == GF2_INVALID_ARGUMENT);
    CHECK(binary_matrix_multiply(&a, &b, &b) == GF2_INVALID_ARGUMENT);
    CHECK(memcmp(c.data, expected, sizeof(expected)) == 0);
    binary_matrix_free(&a);
    binary_matrix_free(&b);
    CHECK(binary_matrix_init(&a, 2, 2) == GF2_OK);
    CHECK(binary_matrix_init(&b, 2, 2) == GF2_OK);
    for (unsigned x = 0; x < 16; ++x) {
        for (unsigned y = 0; y < 16; ++y) {
            for (size_t k = 0; k < 4; ++k) {
                a.data[k] = (uint8_t)((x >> k) & 1u);
                b.data[k] = (uint8_t)((y >> k) & 1u);
            }
            CHECK(binary_matrix_multiply(&a, &b, &c) == GF2_OK);
            for (size_t i = 0; i < 2; ++i) {
                for (size_t j = 0; j < 2; ++j) {
                    unsigned sum = a.data[i * 2] * b.data[j] +
                                   a.data[i * 2 + 1] * b.data[2 + j];
                    CHECK(c.data[i * 2 + j] == sum % 2);
                }
            }
        }
    }
    CHECK(binary_matrix_multiply(&a, &a, &c) == GF2_OK);
    binary_matrix_free(&a);
    binary_matrix_free(&b);
    binary_matrix_free(&c);
}

int main(void)
{
    test_storage();
    test_load();
    test_syndrome();
    test_boundaries();
    test_matrix_product();
    printf("GF(2) checks passed (%u checks; exhaustive 16-vector oracle)\n", checks);
    return EXIT_SUCCESS;
}
