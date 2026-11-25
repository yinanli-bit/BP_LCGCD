#ifndef BPLCGCD_SIMULATION_H
#define BPLCGCD_SIMULATION_H

#include "bplcgcd/decoder.h"

typedef struct {
    double ebn0_db;
    uint64_t seed, max_frames, max_frame_errors; 
} SimulationConfig;

typedef struct {
    size_t n, k;
    double sigma2;
    uint64_t frames, frame_errors, bit_errors, hd_bit_errors;
    uint64_t proven_frames, truncated_frames;
    DecoderStats totals; /* Sums of per-frame counters, not batch/query aliases. */
    bool error_limit_hit; /* Paired runs: common stop because either decoder hit the cap. */
    uint64_t early_stop_counts[DECODE_STOP_REASON_COUNT];
    uint64_t p50_queries, p95_queries, p99_queries; //逐帧完整候选数，nearest-rank 分位数
    double runtime_seconds; //已计入帧的译码调用进程 CPU 时间，包含内部 printf
    bool runtime_available; //clock() 不可用时，CSV=NA，JSON=null
} SimulationResult;

/* Single-thread, one Eb/N0, all-zero classical transmission. */
DecodeStatus simulation_run(const BinaryMatrix *h, const BinarySystematic *form,
                            const DecodeConfig *decoder, const SimulationConfig *config,
                            SimulationResult *result);

/* Single-thread, shared all-zero codeword/y/LLR per frame. lc_config names
 * lcgcd; its candidate_limit applies to BOTH algorithms. results[0]=GCD,
 * results[1]=LC-GCD. Stop together after either error cap or max_frames.
 * Optional borrowed writable paired_trace prints exact hex samples/results.
 * No frame is counted unless both algorithms return a valid estimate. */
DecodeStatus simulation_compare(const BinaryMatrix *h, const BinarySystematic *form,
                               const DecodeConfig *lc_config, const SimulationConfig *config,
                               SimulationResult results[2], FILE *paired_trace);

/* Version 1 schema: header + one completed-run row. BER uses all N codeword bits.
 * Caller owns/closes the writable stream; returns 0 on argument or write error. */
int simulation_write_csv(FILE *stream, const DecodeConfig *decoder,
                         const SimulationConfig *config, const SimulationResult *result);
/* Same versioned columns; header and two paired rows. */
int simulation_compare_write_csv(FILE *stream, const DecodeConfig *lc_config,
                                 const SimulationConfig *config, const SimulationResult results[2]);
/* Same versioned fields as CSV; a JSON array of one or two run records. */
int simulation_write_json(FILE *stream, const DecodeConfig *decoder,
                          const SimulationConfig *config, const SimulationResult *result);
int simulation_compare_write_json(FILE *stream, const DecodeConfig *lc_config,
                                  const SimulationConfig *config, const SimulationResult results[2]);

#endif
