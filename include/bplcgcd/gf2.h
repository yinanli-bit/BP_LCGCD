#ifndef BPLCGCD_GF2_H //当宏 X 未定义时才编译中间代码，保护头文件，避免重复定义
#define BPLCGCD_GF2_H

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

typedef enum {
    GF2_OK = 0,
    GF2_INVALID_ARGUMENT,
    GF2_DIMENSION_MISMATCH,
    GF2_FORMAT_ERROR,
    GF2_IO_ERROR,
    GF2_ALLOCATION_FAILED,
    GF2_INCONSISTENT_SYNDROME
} Gf2Status;


typedef struct {
    size_t rows;
    size_t cols;
    uint8_t *data; /* Entry (i, j) is data[i * cols + j], with zero-based indices. */
} BinaryMatrix;

typedef struct {
    size_t length;
    uint8_t *data; /* Entry j corresponds to matrix column j. */
} BitVector;


Gf2Status binary_matrix_init(BinaryMatrix *matrix, size_t rows, size_t cols);
void binary_matrix_free(BinaryMatrix *matrix);
Gf2Status bit_vector_init(BitVector *vector, size_t length);
void bit_vector_free(BitVector *vector);


Gf2Status binary_matrix_load(FILE *stream, BinaryMatrix *matrix);


Gf2Status binary_matrix_multiply_vector(const BinaryMatrix *matrix,const BitVector *vector,BitVector *result);

Gf2Status binary_matrix_multiply(const BinaryMatrix *left,const BinaryMatrix *right,BinaryMatrix *result);

/* Owning offline result, initialized with {0}; members are read-only to callers.
 * reduced = U * H * Q = [I_rank P; 0 0], retaining all original rows.
 * The top rank rows are H_sys; K = reduced.cols - rank; either may be zero.
 * U = row_transform; Q[sys_to_orig[j], j] = 1. Row XORs require a full U,
 * not just a row permutation. H itself is never modified. */
typedef struct {
    size_t rank;
    BinaryMatrix reduced;
    BinaryMatrix row_transform;
    size_t *sys_to_orig; /* e_sys[j] = e_orig[sys_to_orig[j]]. */
    size_t *orig_to_sys; /* e_orig[j] = e_sys[orig_to_sys[j]]. */
} BinarySystematic;

/* Build requires an empty {0} result; failure leaves it unchanged. */
Gf2Status binary_matrix_systematize(const BinaryMatrix *matrix,BinarySystematic *form);
void binary_systematic_free(BinarySystematic *form);

//交织error[]
Gf2Status binary_systematic_to_system(const BinarySystematic *form,const BitVector *input, BitVector *result);
Gf2Status binary_systematic_to_original(const BinarySystematic *form,const BitVector *input, BitVector *result);

Gf2Status binary_systematic_transform_syndrome(const BinarySystematic *form,const BitVector *input,BitVector *result);

#endif
