#include "bplcgcd/gf2.h"

#include <ctype.h>
#include <stdlib.h>

Gf2Status binary_matrix_init(BinaryMatrix *matrix, size_t rows, size_t cols)
{
    if (!matrix || matrix->data || matrix->rows || matrix->cols ||
        rows == 0 || cols == 0 || rows > SIZE_MAX / cols) {
        return GF2_INVALID_ARGUMENT;
    }
    uint8_t *data = calloc(rows * cols, sizeof(*data));
    if (!data) {
        return GF2_ALLOCATION_FAILED;
    }
    *matrix = (BinaryMatrix){rows, cols, data};
    return GF2_OK;
}

void binary_matrix_free(BinaryMatrix *matrix)
{
    if (matrix) {
        free(matrix->data);
        *matrix = (BinaryMatrix){0};
    }
}

Gf2Status bit_vector_init(BitVector *vector, size_t length)
{
    if (!vector || vector->data || vector->length || length == 0) {
        return GF2_INVALID_ARGUMENT;
    }
    uint8_t *data = calloc(length, sizeof(*data));
    if (!data) {
        return GF2_ALLOCATION_FAILED;
    }
    *vector = (BitVector){length, data};
    return GF2_OK;
}

void bit_vector_free(BitVector *vector)
{
    if (vector) {
        free(vector->data);
        *vector = (BitVector){0};
    }
}

static int next_nonspace(FILE *stream)
{
    int ch;
    do {
        ch = fgetc(stream);
    } while (ch != EOF && isspace((unsigned char)ch));
    return ch;
}

static Gf2Status read_dimension(FILE *stream, size_t *dimension)
{
    int ch = next_nonspace(stream);
    size_t value = 0;
    if (ch < '0' || ch > '9') {
        return ferror(stream) ? GF2_IO_ERROR : GF2_FORMAT_ERROR;
    }
    do {
        size_t digit = (size_t)(ch - '0');
        if (value > (SIZE_MAX - digit) / 10) {
            return GF2_FORMAT_ERROR;
        }
        value = value * 10 + digit;//可能是多位数
        ch = fgetc(stream);
    } while (ch >= '0' && ch <= '9');
    if (ferror(stream)) {
        return GF2_IO_ERROR;
    }
    if (value == 0 || (ch != EOF && !isspace((unsigned char)ch))) {
        return GF2_FORMAT_ERROR;
    }
    *dimension = value;
    return GF2_OK;
}

// 矩阵文件格式，首行两个数字 行数 和 列数
Gf2Status binary_matrix_load(FILE *stream, BinaryMatrix *matrix)
{
    if (!stream || !matrix || matrix->data || matrix->rows || matrix->cols) {
        return GF2_INVALID_ARGUMENT;
    }
    size_t rows, cols;
    Gf2Status status = read_dimension(stream, &rows);
    if (status != GF2_OK) {
        return status;
    }
    status = read_dimension(stream, &cols);
    if (status != GF2_OK) {
        return status;
    }
    if (rows > SIZE_MAX / cols) {
        return GF2_FORMAT_ERROR;
    }
    BinaryMatrix loaded = {0};
    status = binary_matrix_init(&loaded, rows, cols);
    if (status != GF2_OK) {
        return status;
    }
    for (size_t i = 0; i < rows * cols; ++i) {
        int ch = next_nonspace(stream);
        if (ch != '0' && ch != '1') {
            status = GF2_FORMAT_ERROR;
            break;
        }
        loaded.data[i] = (uint8_t)(ch - '0');
        ch = fgetc(stream);
        if (ch != EOF && !isspace((unsigned char)ch)) {
            status = GF2_FORMAT_ERROR;
            break;
        }
    }
    if (status == GF2_OK && next_nonspace(stream) != EOF) {
        status = GF2_FORMAT_ERROR;
    }
    if (ferror(stream)) {
        status = GF2_IO_ERROR;
    }
    if (status != GF2_OK) {
        binary_matrix_free(&loaded);
        return status;
    }
    *matrix = loaded;
    return GF2_OK;
}

Gf2Status binary_matrix_multiply_vector(const BinaryMatrix *matrix,const BitVector *vector,BitVector *result)
{
    if (!matrix || !vector || !result || !matrix->data || !vector->data ||
        !result->data || matrix->rows == 0 || matrix->cols == 0 ||
        matrix->rows > SIZE_MAX / matrix->cols ||
        vector->length == 0 || result->length == 0 ||
        result->data == vector->data || result->data == matrix->data) {
        return GF2_INVALID_ARGUMENT;
    }
    if (vector->length != matrix->cols || result->length != matrix->rows) {
        return GF2_DIMENSION_MISMATCH;
    }
    /* Validate every input before writing any result, including late errors. */
    for (size_t i = 0; i < matrix->rows * matrix->cols; ++i) {
        if (matrix->data[i] > 1) {
            return GF2_INVALID_ARGUMENT;
        }
    }
    for (size_t j = 0; j < vector->length; ++j) {
        if (vector->data[j] > 1) {
            return GF2_INVALID_ARGUMENT;
        }
    }
    for (size_t i = 0; i < matrix->rows; ++i) {
        uint8_t parity = 0;
        for (size_t j = 0; j < matrix->cols; ++j) {
            parity ^= matrix->data[i * matrix->cols + j] & vector->data[j];
        }
        result->data[i] = parity;
    }
    return GF2_OK;
}

Gf2Status binary_matrix_multiply(const BinaryMatrix *left,const BinaryMatrix *right,BinaryMatrix *result)
{
    if (!left || !right || !result || !left->data || !right->data ||
        !result->data || !left->rows || !left->cols || !right->rows ||
        !right->cols || !result->rows || !result->cols ||
        left->rows > SIZE_MAX / left->cols ||
        right->rows > SIZE_MAX / right->cols ||
        result->rows > SIZE_MAX / result->cols ||
        result->data == left->data || result->data == right->data) {
        return GF2_INVALID_ARGUMENT;
    }
    if (left->cols != right->rows || result->rows != left->rows ||
        result->cols != right->cols) {
        return GF2_DIMENSION_MISMATCH;
    }
    for (size_t i = 0; i < left->rows * left->cols; ++i) {
        if (left->data[i] > 1) return GF2_INVALID_ARGUMENT;
    }
    for (size_t i = 0; i < right->rows * right->cols; ++i) {
        if (right->data[i] > 1) return GF2_INVALID_ARGUMENT;
    }
    for (size_t i = 0; i < result->rows; ++i) {
        for (size_t j = 0; j < result->cols; ++j) {
            uint8_t parity = 0;
            for (size_t k = 0; k < left->cols; ++k) {
                parity ^= left->data[i * left->cols + k] &
                          right->data[k * right->cols + j];
            }
            result->data[i * result->cols + j] = parity;
        }
    }
    return GF2_OK;
}
