#ifndef BPLCGCD_QUANTUM_CHANNEL_H
#define BPLCGCD_QUANTUM_CHANNEL_H

#include "bplcgcd/quantum.h"

typedef struct {
    double p_i, p_x, p_y, p_z;
} PauliProbabilities;

QuantumStatus quantum_depolarizing_probabilities(double p, PauliProbabilities *prob);
/* Strictly positive normalized probabilities -> finite natural prior XYZ LLR.
 * Caller may invoke per qubit to support nonuniform/biased Pauli channels. */
QuantumStatus quantum_prior_llr(const PauliProbabilities *prob, double prior_xyz[3]);
/* One uniform draw per qubit. Stateless frame replay; no allocation.
 * Domain-separated from classical AWGN. Sampling supports p=0 and p=1. */
QuantumStatus quantum_depolarizing_sample(double p, uint64_t seed,
                                          uint64_t frame_index,
                                          size_t length, Pauli *error);

#endif
