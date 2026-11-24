#ifndef BPLCGCD_LOCAL_CONSTRAINTS_H
#define BPLCGCD_LOCAL_CONSTRAINTS_H

#include "bplcgcd/gf2.h"

/* Caller-owned buffers. basis_order has m entries in ORIGINAL column indices:
 * first delta are A, the rest are Ac, ordered by decreasing flip weight, then
 * increasing original index. P-block column j maps to sys_to_orig[m+j]. */
typedef struct {
    size_t *basis_order;
    BinaryMatrix p_a, p_ac; /* delta x K and (m-delta) x K. */
    BitVector s_a, s_ac;   /* delta and m-delta entries of U*s_res. */
} LocalConstraints;


//实现论文公式： f_I[A]  XOR P_A · f_Pᵀ  = s_A
//              f_I[Ac] XOR P_Ac · f_Pᵀ = s_Ac
Gf2Status local_constraints_prepare(const BinarySystematic *form,
                                    const double *flip_weights, size_t length,
                                    const BitVector *s_res, size_t delta,
                                    LocalConstraints *result, BitVector *syndrome_scratch);

#endif
