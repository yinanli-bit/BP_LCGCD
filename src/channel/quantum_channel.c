#include "bplcgcd/quantum_channel.h"

#include <math.h>

QuantumStatus quantum_depolarizing_probabilities(double p, PauliProbabilities *prob)
{
    if (!prob || !isfinite(p) || p < 0 || p > 1) 
        return QUANTUM_INVALID_ARGUMENT;
    *prob = (PauliProbabilities){1-p, p/3, p/3, p/3};
    return QUANTUM_OK;
}

QuantumStatus quantum_prior_llr(const PauliProbabilities *prob, double prior_xyz[3])
{
    if (!prob || !prior_xyz) 
        return QUANTUM_INVALID_ARGUMENT;
    const double p[] = {prob->p_i, prob->p_x, prob->p_y, prob->p_z};
    double sum = 0;
    for (size_t k = 0; k < 4; ++k) {
        if (!isfinite(p[k]) || p[k] <= 0 || p[k] > 1) 
            return QUANTUM_INVALID_ARGUMENT;
        sum += p[k];
    }
    if (fabs(sum - 1) > 1e-12) 
        return QUANTUM_INVALID_ARGUMENT;
    /* Difference of logs also works when the probability ratio overflows. */
    for (size_t k = 0; k < 3; ++k) 
        prior_xyz[k] = log(p[0]) - log(p[k+1]);
    return QUANTUM_OK;
}

/* Same fixed SplitMix64 mixer as the existing classical sampler. */
static uint64_t mix64(uint64_t value)
{
    value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31);
}

//作为模拟的真实错误，与后面译码结果进行对照
QuantumStatus quantum_depolarizing_sample(double p, uint64_t seed,
                                          uint64_t frame_index,
                                          size_t length, Pauli *error)
{
    if (!error || !length || !isfinite(p) || p < 0 || p > 1)
        return QUANTUM_INVALID_ARGUMENT;
    uint64_t state = seed ^ mix64(frame_index) ^ UINT64_C(0x5155414e54554d31);
    for (size_t c = 0; c < length; ++c) {
        state += UINT64_C(0x9e3779b97f4a7c15);
        double u = ((double)(mix64(state) >> 12) + 0.5) * 0x1p-52;
        if (u < 1-p) error[c] = PAULI_I;
        else if (u < 1-p + p/3) error[c] = PAULI_X;
        else if (u < 1-p + 2*(p/3)) error[c] = PAULI_Y;
        else error[c] = PAULI_Z;
    }
    return QUANTUM_OK;
}
