#include "bplcgcd/local_constraints.h"

#include <math.h>
#include <string.h>

void print_local_constraints(const LocalConstraints *local)
{
    printf("basis_order: ");
    for (size_t i = 0; i < local->p_a.rows + local->p_ac.rows; ++i) {
        printf("%zu ", local->basis_order[i]);
    }
    printf("\n");
    printf("p_a:\n");
    for (size_t i = 0; i < local->p_a.rows; ++i) {
        for (size_t j = 0; j < local->p_a.cols; ++j) {
            printf("%d ", local->p_a.data[i * local->p_a.cols + j]);
        }
        printf("\n");
    }
    printf("p_ac:\n");
    for (size_t i = 0; i < local->p_ac.rows; ++i) {
        for (size_t j = 0; j < local->p_ac.cols; ++j) {
            printf("%d ", local->p_ac.data[i * local->p_ac.cols + j]);
        }
        printf("\n");
    }
    printf("s_a: ");
    for (size_t i = 0; i < local->s_a.length; ++i) {
        printf("%d ", local->s_a.data[i]);
    }
    printf("\n");
    printf("s_ac: ");
    for (size_t i = 0; i < local->s_ac.length; ++i) {
        printf("%d ", local->s_ac.data[i]);
    }
    printf("\n");
}

Gf2Status local_constraints_prepare(const BinarySystematic *form,
                                    const double *flip_weights, size_t length,
                                    const BitVector *s_res, size_t delta,
                                    LocalConstraints *result, BitVector *syndrome_scratch)
{
    if (!form || !form->reduced.data || !flip_weights || !s_res || !result || !syndrome_scratch) {
        return GF2_INVALID_ARGUMENT;
    }
    size_t n = form->reduced.cols, m = form->rank, rows = form->reduced.rows;
    if (m > n || m > rows || length != n || delta > m) return GF2_DIMENSION_MISMATCH;

    size_t k = n - m, rest = m - delta;

    if (s_res->length != rows || syndrome_scratch->length != rows ||result->p_a.rows != delta || result->p_a.cols != k ||result->p_ac.rows != rest || result->p_ac.cols != k ||result->s_a.length != delta || result->s_ac.length != rest) 
            return GF2_DIMENSION_MISMATCH;
    if ((m && !result->basis_order) || (delta && !result->s_a.data) || (rest && !result->s_ac.data) ||(delta && k && !result->p_a.data) || (rest && k && !result->p_ac.data)) 
            return GF2_INVALID_ARGUMENT;
    for (size_t j = 0; j < n; ++j) {
        if (!isfinite(flip_weights[j]) || flip_weights[j] < 0) return GF2_INVALID_ARGUMENT;
    }
    Gf2Status status = binary_systematic_transform_syndrome(form, s_res, syndrome_scratch);
    if (status != GF2_OK) return status;

    /* Only the fixed basis participates, even if a free column is more reliable. */
    //完成对于basis的排序，对应original的列索引
    for (size_t i = 0; i < m; ++i) {
        size_t original = form->sys_to_orig[i], pos = i;
        while (pos) {
            size_t previous = result->basis_order[pos - 1];
            if (flip_weights[previous] > flip_weights[original] ||(flip_weights[previous] == flip_weights[original] && previous < original)) 
                break;
            result->basis_order[pos] = previous;
            --pos;
        }
        result->basis_order[pos] = original;
    } 
    for (size_t part = 0; part < 2; ++part) {
        BinaryMatrix *p = part ? &result->p_ac : &result->p_a;
        BitVector *s = part ? &result->s_ac : &result->s_a;
        size_t start = part ? delta : 0; //可靠位是前delta位，剩余位是后m-delta位
        for (size_t i = 0; i < s->length; ++i) {
            size_t basis = form->orig_to_sys[result->basis_order[start + i]]; //索引按系统形式顺序
            if (k) memcpy(p->data + i * k, form->reduced.data + basis * n + m, k); //将系统化矩阵的P块复制到p_a或p_ac中
            s->data[i] = syndrome_scratch->data[basis];
        }
    }
    print_local_constraints(result);
    return GF2_OK;
}
