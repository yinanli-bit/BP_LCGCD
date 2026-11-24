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
    bool error_limit_hit;
} SimulationResult;

/* Single-thread, one Eb/N0, all-zero classical transmission. */
DecodeStatus simulation_run(const BinaryMatrix *h, const BinarySystematic *form,
                            const DecodeConfig *decoder, const SimulationConfig *config,
                            SimulationResult *result);

/* Writes a header + one completed-run row. BER uses all N codeword bits.
 * Caller owns/closes the writable stream; returns 0 on argument or write error. */
int simulation_write_csv(FILE *stream, const DecodeConfig *decoder,
                         const SimulationConfig *config, const SimulationResult *result);

#endif
