#ifndef BPLCGCD_QUANTUM_DECODER_H
#define BPLCGCD_QUANTUM_DECODER_H

#include "bplcgcd/decoder.h"
#include "bplcgcd/qbp.h"
#include "bplcgcd/slva.h"

typedef struct { double w[4]; } PauliSwt; /* Relative I/X/Y/Z weights; min=0. */
typedef struct { double x[2], z[2]; } PartialSwt;

/* Caller-owned arrays of posterior->length entries, natural logarithms.
 * These helpers transform POSTERIOR information, never channel parameters. */
QuantumStatus qbp_posterior_to_swt(const QbpResult *posterior, PauliSwt *swt);
/* hard/reliability order is (EX|EZ), length 2*n. Partial weights may be negative
 * because SWT has a relative normalization. Output buffers must not overlap. */
QuantumStatus quantum_swt_to_partial(const PauliSwt *swt, size_t n,
                                     PartialSwt *partial, BitVector *hard,
                                     double *reliability);

typedef enum { BP_POST_ON_FAILURE = 0, BP_POST_ALWAYS } BpPostPolicy;

typedef enum {
    QUANTUM_STOP_TRIVIAL = 0, /* Default: rigorous Pauli lower bound, Eq. (23). */
    QUANTUM_STOP_DISABLED,
    QUANTUM_STOP_DAI         /* Eq. (24)-(26); approximate stop is not a proof. */
} QuantumStopPolicy;

typedef enum {
    QUANTUM_EARLY_NONE = 0, QUANTUM_EARLY_TRIVIAL, QUANTUM_EARLY_DAI
} QuantumEarlyStop;

typedef struct {
    size_t delta;
    uint64_t candidate_limit; /* Complete Pauli candidates; 0=unlimited. */
    size_t state_limit, slva_max_candidates, slva_max_bytes; /* 0: classical defaults. */
    BpPostPolicy post_policy; /* Default zero: BP_POST_ON_FAILURE. */
    QuantumStopPolicy stop_policy;
} QuantumLcgcdConfig;

/* Initialize once with {0}. Owning workspace; members are read-only to callers.
 * code is borrowed and must outlive this object; config is copied.
 * Hsyn and S are systematized only here, never in frame/search routines. */
typedef struct {
    const QuantumCode *code;
    QuantumLcgcdConfig config;
    BinarySystematic syndrome_form, stabilizer_form;
    QbpWorkspace bp;
    PauliSwt *swt;
    PartialSwt *partial_swt;
    double *reliability;
    Pauli *candidate, *recovery;
    BitVector reference, residual, transformed, system, partial, binary, syndrome, local_rhs;
    LocalConstraints local;
    uint8_t *local_p;
    SyndromeTrellis trellis;
    SlvaGenerator list;
    long double stop_correction; /* Pauli coupling + conservative sum roundoff. */
    double dai_threshold;
} QuantumDecoderWorkspace;

/* Recovery is borrowed until the next decode/free. Standalone posterior views
 * retain the caller's input buffers; BP-wrapper views belong to workspace.bp.
 * A valid result always matches measured syndrome; logical success is separate.
 * optimality_proven means a Pauli bound or exhaustion proves the minimum of
 * the supplied floating posterior score, never coset/physical ML. DAI-only
 * termination has proof=false and early_stop=QUANTUM_EARLY_DAI. */
typedef struct {
    DecodeStatus status;
    size_t length;
    const Pauli *recovery;
    double cost;
    bool candidate_valid, syndrome_matched, optimality_proven, postprocess_invoked;
    QbpResult posterior;
    DecoderStats stats;
    DecodeStopReason stop_reason;
    QuantumEarlyStop early_stop; /* DAI uses stop_reason=NONE, not a proven bound. */
    double remaining_lower_bound, dai_threshold;
} QuantumDecodeResult;

typedef struct { 
    bool syndrome_matched, logical_success; 
} QuantumRecoveryCheck;

/* Resource-cap failure returns DECODE_TRUNCATED and leaves workspace empty. */
DecodeStatus quantum_decoder_init(QuantumDecoderWorkspace *workspace,
                                  const QuantumCode *code, const QuantumLcgcdConfig *config);
void quantum_decoder_free(QuantumDecoderWorkspace *workspace);
/* Standalone post-processing of supplied posterior, including failed BP output.
 * Caller-owned posterior buffers must outlive the returned posterior view. */
DecodeStatus quantum_lcgcd_decode(QuantumDecoderWorkspace *workspace,
                                  const QbpResult *posterior, const BitVector *measured,
                                  QuantumDecodeResult *result);
/* Single QNMS run; post-process according to workspace.config.post_policy. */
DecodeStatus quantum_bp_lcgcd_decode(QuantumDecoderWorkspace *workspace,
                                     const QbpConfig *bp_config, const QbpPrior *prior,
                                     const BitVector *measured, QuantumDecodeResult *result);
/* E/R in original qubit order. Uses the offline S row space and reusable scratch;
 * no allocation/GE. Recovery need not match syndrome; both flags are reported. */
QuantumStatus quantum_recovery_check(QuantumDecoderWorkspace *workspace,
                                     const Pauli *physical_error, const Pauli *recovery,
                                     size_t length, QuantumRecoveryCheck *result);

#endif
