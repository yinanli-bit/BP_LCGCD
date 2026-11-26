#ifndef BPLCGCD_QBP_H
#define BPLCGCD_QBP_H

#include "bplcgcd/quantum.h"

typedef struct {
    size_t length;
    const double *prior_x, *prior_y, *prior_z;
} QbpPrior;

typedef struct {
    size_t max_iterations; /* Must be positive. */
    double nms_scale;      /* Direct multiplication, 0 < scale <= 1. */
    double message_clip;   /* Finite positive C2V magnitude bound before scale;
                            * zero selects 50. Numerical guard, not quantization. */
} QbpConfig;

/* Owning, initialized with {0}; code must outlive workspace. Do not copy or
 * mutate members. Independent workspaces can share the same read-only code. */
typedef struct {
    const QuantumCode *code;
    double *messages;
    double *v2c_x, *v2c_y, *v2c_z, *c2v, *commutation;
    double *post_x, *post_y, *post_z;
    Pauli *hard_estimate;
    uint8_t *syndrome;
} QbpWorkspace;

/* Borrowed output, valid until next decode/free of this workspace. */
typedef struct {
    size_t length, iterations_used;
    int syndrome_matched; /* Convergence only, not logical success. */
    const Pauli *hard_estimate;
    const double *post_x, *post_y, *post_z;
} QbpResult;

QuantumStatus qbp_workspace_init(QbpWorkspace *workspace, const QuantumCode *code);
void qbp_workspace_free(QbpWorkspace *workspace);
/* 2025 Eq. (6); msg_xyz must contain finite X/Y/Z messages, stabilizer must
 * be X/Y/Z. Invalid arguments return NAN. Stable log-sum-exp evaluation. */
double qbp_commutation_llr(const double msg_xyz[3], Pauli stabilizer_pauli);
/* Floating flooding QNMS, 2025 Eq. (3)-(10), (12). No allocation or topology
 * mutation. Prior buffers must not overlap workspace storage. QUANTUM_OK
 * includes iteration exhaustion; inspect result.syndrome_matched separately.
 * Error clears result; workspace is reusable after errors. */
QuantumStatus qbp_decode(QbpWorkspace *workspace, const QbpConfig *config,
                         const QbpPrior *prior, const BitVector *measured_syndrome,
                         QbpResult *result);

#endif
