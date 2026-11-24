#include "bplcgcd/simulation.h"
#include "bplcgcd/channel.h"

#include <inttypes.h>
#include <stdlib.h>
#include <string.h>

DecodeStatus simulation_run(const BinaryMatrix *h, const BinarySystematic *form,
                            const DecodeConfig *decoder, const SimulationConfig *config,
                            SimulationResult *result)
{
    if (!h || !h->data || !h->cols || !form || form->rank >= h->cols || !decoder ||
        !decoder->decoder_name || !config || !result || !config->max_frames ||
        !config->max_frame_errors || h->cols > SIZE_MAX / sizeof(double) ||
        config->max_frames > UINT64_MAX / h->cols) return DECODE_INVALID_ARGUMENT;
    if (strcmp(decoder->decoder_name, "gcd")) return DECODE_NOT_IMPLEMENTED;
    if (form->reduced.rows != h->rows || form->reduced.cols != h->cols) return DECODE_DIMENSION_MISMATCH;
    size_t n = h->cols, k = n - form->rank;
    double sigma2;
    ChannelStatus channel = bpsk_noise_variance((double)k / n, config->ebn0_db, &sigma2);
    if (channel != CHANNEL_OK) return channel == CHANNEL_NUMERIC_ERROR ? DECODE_NUMERIC_ERROR : DECODE_INVALID_ARGUMENT;
    *result = (SimulationResult){.n=n, .k=k, .sigma2=sigma2};
    BitVector zero = {0};
    DecodeResult decoded = {0};
    double *received = malloc(n * sizeof(*received)), *llr = malloc(n * sizeof(*llr));
    DecodeStatus status = DECODE_ALLOCATION_FAILED;
    if (!received || !llr || bit_vector_init(&zero, n) != GF2_OK ||
        bit_vector_init(&decoded.estimate, n) != GF2_OK) goto finish;
    DecodeInput input = {.mode=DECODE_CLASSICAL, .h=h, .length=n, .llr=llr, .systematic=form};
    while (result->frames < config->max_frames && result->frame_errors < config->max_frame_errors) {
        channel = bpsk_awgn_transmit(&zero, sigma2, config->seed, result->frames, received);
        if (channel == CHANNEL_OK) 
            channel = bpsk_compute_llr(received, n, sigma2, llr);
        if (channel != CHANNEL_OK) 
            { status = DECODE_NUMERIC_ERROR; goto finish; }
        status = decoder_decode(decoder, &input, &decoded);
        if (status != DECODE_SUCCESS && status != DECODE_TRUNCATED) goto finish;
        if (!decoded.candidate_valid) { status = DECODE_INVALID_CANDIDATE; goto finish; }
        uint64_t errors = 0;
        for (size_t j = 0; j < n; ++j) {
            errors += decoded.estimate.data[j]; /* Compare against transmitted zero. */
            result->hd_bit_errors += llr[j] < 0;
        }
        ++result->frames;
        result->bit_errors += errors;
        result->frame_errors += errors != 0;
        result->truncated_frames += status == DECODE_TRUNCATED;
        result->proven_frames += decoded.optimality_proven;
        result->totals.candidates_popped += decoded.stats.candidates_popped;
        result->totals.candidates_evaluated += decoded.stats.candidates_evaluated;
        result->totals.nodes_expanded += decoded.stats.nodes_expanded;
        result->totals.batches += decoded.stats.batches;
    }
    result->error_limit_hit = result->frame_errors >= config->max_frame_errors;
    status = DECODE_SUCCESS; /* The run completed; proof/truncation remain per frame. */
finish:
    bit_vector_free(&decoded.estimate);
    bit_vector_free(&zero);
    free(llr);
    free(received);
    return status;
}

int simulation_write_csv(FILE *stream, const DecodeConfig *decoder,
                         const SimulationConfig *config, const SimulationResult *result)
{
    if (!stream || !decoder || !config || !result || !result->frames) return 0;
    fputs("ebn0_db,seed,n,k,rate,sigma2,max_frames,max_frame_errors,candidate_limit,"
          "frames,frame_errors,bit_errors,hd_bit_errors,codeword_ber,fer,hd_codeword_ber,"
          "candidates_popped,candidates_evaluated,nodes_expanded,batches,avg_candidates_evaluated,"
          "proven_frames,truncated_frames,stop_reason\n", stream);
    double bits = (double)result->frames * result->n;
    int written = fprintf(stream, "%.17g,%" PRIu64 ",%zu,%zu,%.17g,%.17g,%" PRIu64 ",%" PRIu64 ",%" PRIu64
        ",%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%.17g,%.17g,%.17g"
        ",%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%.17g,%" PRIu64 ",%" PRIu64 ",%s\n",
        config->ebn0_db, config->seed, result->n, result->k, (double)result->k / result->n,
        result->sigma2, config->max_frames, config->max_frame_errors, decoder->candidate_limit,
        result->frames, result->frame_errors, result->bit_errors, result->hd_bit_errors,
        result->bit_errors / bits, (double)result->frame_errors / result->frames, result->hd_bit_errors / bits,
        result->totals.candidates_popped, result->totals.candidates_evaluated,
        result->totals.nodes_expanded, result->totals.batches,
        (double)result->totals.candidates_evaluated / result->frames,
        result->proven_frames, result->truncated_frames, result->error_limit_hit ? "error_limit" : "frame_limit");
    return written >= 0 && !ferror(stream);
}
