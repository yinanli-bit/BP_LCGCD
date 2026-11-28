#include "bplcgcd/quantum_decoder.h"

#include <math.h>

QuantumStatus qbp_posterior_to_swt(const QbpResult *post, PauliSwt *swt)
{
    if (!post || !post->length || !post->post_x || !post->post_y || !post->post_z || !swt)
        return QUANTUM_INVALID_ARGUMENT;
    for (size_t i = 0; i < post->length; ++i) {
        const double relative[] = {0,post->post_x[i],post->post_y[i],post->post_z[i]};
        double minimum = 0;
        for (size_t p = 1; p < 4; ++p) {
            if (!isfinite(relative[p])) return QUANTUM_INVALID_ARGUMENT;
            minimum = fmin(minimum, relative[p]);
        }
        /* 2025 Eq. (21): log(max q)-log(qP). Constants cancel in full scoring. */
        for (size_t p = 0; p < 4; ++p) {
            swt[i].w[p] = relative[p]-minimum;
            if (!isfinite(swt[i].w[p])) return QUANTUM_NUMERIC_ERROR;
        }
    }
    return QUANTUM_OK;
}

static double soft_add(double a, double b)
{
    double low = fmin(a,b), high = fmax(a,b);
    return low-log1p(exp(low-high)); /* -log(exp(-a)+exp(-b)), stable. */
}

QuantumStatus quantum_swt_to_partial(const PauliSwt *swt, size_t n,
                                     PartialSwt *partial, BitVector *hard,
                                     double *reliability)
{
    if (!swt || !n || n > SIZE_MAX/2 || !partial || !hard || !hard->data || !reliability)
        return QUANTUM_INVALID_ARGUMENT;
    if (hard->length != 2*n) 
        return QUANTUM_DIMENSION_MISMATCH;
    for (size_t i = 0; i < n; ++i) {
        const double *w = swt[i].w;
        for (size_t p = 0; p < 4; ++p)
            if (!isfinite(w[p]) || w[p] < 0) return QUANTUM_INVALID_ARGUMENT;
        /* 2024 Eq. (13)-(15): marginalize Pauli masses, including Y in both 1s. */
        partial[i] = (PartialSwt){
            .x={soft_add(w[PAULI_I],w[PAULI_Z]), soft_add(w[PAULI_X],w[PAULI_Y])},
            .z={soft_add(w[PAULI_I],w[PAULI_X]), soft_add(w[PAULI_Z],w[PAULI_Y])}};
        hard->data[i] = partial[i].x[1] < partial[i].x[0];
        hard->data[n+i] = partial[i].z[1] < partial[i].z[0];
        reliability[i] = fabs(partial[i].x[1]-partial[i].x[0]);
        reliability[n+i] = fabs(partial[i].z[1]-partial[i].z[0]);
    }
    return QUANTUM_OK;
}
