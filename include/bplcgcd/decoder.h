#ifndef BPLCGCD_DECODER_H
#define BPLCGCD_DECODER_H

#include "bplcgcd/gf2.h"
#include <stdbool.h>

typedef enum { DECODE_CLASSICAL, DECODE_SYNDROME, DECODE_BP_SYNDROME } DecodeMode;
typedef enum {
    DECODE_NOT_RUN = 0,
    DECODE_SUCCESS, /* Operation succeeded; inspect candidate_valid for a result. */
    DECODE_TRUNCATED, /* Budget exhausted; candidate_valid may be false. */
    DECODE_NOT_IMPLEMENTED,
    DECODE_INVALID_ARGUMENT,
    DECODE_DIMENSION_MISMATCH,
    DECODE_INVALID_CANDIDATE,
    DECODE_NUMERIC_ERROR,
    DECODE_ALLOCATION_FAILED,
    DECODE_INCONSISTENT_SYNDROME
} DecodeStatus;

typedef enum {
    DECODE_STOP_NONE = 0,
    DECODE_STOP_BOUND, //早停：未访问候选的部分代价下界不小于当前最佳完整代价
    DECODE_STOP_EXHAUSTED, //搜索空间耗尽，当前候选代价等于最优代价
    DECODE_STOP_BUDGET //搜索空间未耗尽，但候选数达到预算限制
} DecodeStopReason;

/*
 * CLASSICAL: channel llr only; weights, s_target, z_hint must be NULL.
 * SYNDROME: nonnegative weights, s_target, optional z_hint (default zero);
 * weights are flip costs RELATIVE TO z_hint, not costs of final error bits.
 * BP_SYNDROME: signed BP posterior llr and s_target; weights/z_hint are NULL.
 * Both LLR modes derive reference[j]=[llr[j]<0], flip_weights[j]=abs(llr[j]).
 */
typedef struct {
    DecodeMode mode;
    const BinaryMatrix *h;
    size_t length; //length=N=H->cols
    const double *llr;
    const double *weights;
    const BitVector *s_target;
    const BitVector *z_hint;
    /* Required by GCD: borrowed offline S03 form built from this same H.
     * H/form remain unchanged; no per-frame Gaussian elimination. */
    const BinarySystematic *systematic;
} DecodeInput;

typedef struct {
    const char *decoder_name;
    uint64_t candidate_limit; /* Complete candidate evaluations; 0=unlimited. */
    FILE *trace; /* Optional borrowed writable stream; NULL disables GCD trace. */
} DecodeConfig;


//译码过程中的四种计数
typedef struct {
    uint64_t candidates_popped; //弹出数
    uint64_t candidates_evaluated; //完成完整代价计算的候选数
    uint64_t nodes_expanded; //出堆后生成至少一个孩子的父节点数
    uint64_t batches;
} DecoderStats;

/*Read estimate/cost only when candidate_valid. optimality_proven requires a
 * future decoder's proof; check_candidate always leaves it false. */
typedef struct {
    DecodeStatus status;
    BitVector estimate; //is v in CLASSICAL, e_hat otherwise
    double cost; //W(f)=sum(f[j]*flip_weights[j]), relative to the reference.
    bool candidate_valid;
    bool optimality_proven;
    DecoderStats stats;
    DecodeStopReason stop_reason;
} DecodeResult;

/*
 * CLASSICAL: reference=z, s_res=H*z^T. Other modes: reference=z_hint (explicit
 * or derived from BP llr), s_res=s_target XOR H*reference^T.
 * The future common search kernel receives s_res and nonnegative flip_weights,
 * finds H*f^T=s_res and minimizes W(f), without inspecting the mode/reference.
 * S03 checks redundant-row consistency separately. Do not use outputs on error. */
DecodeStatus decode_prepare(const DecodeInput *input, BitVector *reference,BitVector *s_res, double *flip_weights);

DecodeStatus decode_check_candidate(const DecodeInput *input,const double *flip_weights, const BitVector *f,DecodeResult *result, BitVector *syndrome_scratch);
void decode_result_reset(DecodeResult *result); /* Retains buffer, clears metadata. */
void decode_result_mark_truncated(DecodeResult *result); /* Retains best candidate. */
bool decoder_is_known(const char *name);

DecodeStatus decoder_decode(const DecodeConfig *config, const DecodeInput *input,DecodeResult *result);

#endif
