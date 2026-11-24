#include "bplcgcd/simulation.h"
#include "bplcgcd/channel.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

static unsigned checks;
#define CHECK(condition) do { \
    if (!(condition)) { \
        fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, #condition); \
        exit(EXIT_FAILURE); \
    } \
    ++checks; \
} while (0)

/* Independent codeword/distance oracle. For budget=1, the completed empty
 * partial flip keeps the four free hard-decision bits at columns 3..6. */
static void check_oracle(const BinaryMatrix *h, const SimulationConfig *config,
                         const SimulationResult *actual, int budget_one)
{
    uint8_t data[7] = {0};
    BitVector zero = {7, data};
    double y[7], sigma2 = 1 / (2 * (4.0 / 7) * pow(10, config->ebn0_db / 10));
    uint64_t frames = 0, errors = 0, bit_errors = 0, hd_errors = 0;
    while (frames < config->max_frames && errors < config->max_frame_errors) {
        CHECK(bpsk_awgn_transmit(&zero, sigma2, config->seed, frames, y) == CHANNEL_OK);
        double best = INFINITY;
        unsigned chosen = 0;
        for (unsigned word = 0; word < 128; ++word) {
            unsigned invalid = 0;
            for (size_t i = 0; i < h->rows; ++i) {
                unsigned sum = 0;
                for (size_t j = 0; j < 7; ++j) sum += h->data[i * 7 + j] * ((word >> j) & 1u);
                invalid |= sum % 2;
            }
            if (budget_one) {
                for (size_t j = 3; j < 7; ++j) invalid |= ((word >> j) & 1u) != (unsigned)(y[j] < 0);
            }
            if (invalid) continue;
            double distance = 0;
            for (size_t j = 0; j < 7; ++j) {
                double delta = y[j] - (1 - 2 * (int)((word >> j) & 1u));
                distance += delta * delta;
            }
            if (distance < best) { best = distance; chosen = word; }
        }
        CHECK(isfinite(best));
        errors += chosen != 0;
        for (size_t j = 0; j < 7; ++j) {
            bit_errors += (chosen >> j) & 1u;
            hd_errors += y[j] < 0;
        }
        ++frames;
    }
    CHECK(actual->frames == frames && actual->frame_errors == errors);
    CHECK(actual->bit_errors == bit_errors && actual->hd_bit_errors == hd_errors);
    CHECK(actual->error_limit_hit == (errors >= config->max_frame_errors));
    CHECK(actual->n == 7 && actual->k == 4 && actual->sigma2 == sigma2);
}

static void csv_text(const DecodeConfig *decoder, const SimulationConfig *config,
                     const SimulationResult *result, char *buffer, size_t size)
{
    FILE *stream = tmpfile();
    CHECK(stream != NULL);
    CHECK(simulation_write_csv(stream, decoder, config, result));
    CHECK(fflush(stream) == 0);
    rewind(stream);
    size_t count = fread(buffer, 1, size - 1, stream);
    CHECK(!ferror(stream) && feof(stream));
    buffer[count] = '\0';
    CHECK(fclose(stream) == 0);
}

int main(void)
{
    /* Self-authored Hamming [7,4], fourth row redundant: R uses rank, not rows. */
    uint8_t entries[] = {1,0,0,1,1,0,1, 0,1,0,1,0,1,1,
                         0,0,1,0,1,1,1, 1,1,0,0,1,1,0};
    BinaryMatrix h = {4, 7, entries};
    BinarySystematic form = {0};
    CHECK(binary_matrix_systematize(&h, &form) == GF2_OK && form.rank == 3);
    DecodeConfig decoder = {.decoder_name="gcd"};
    SimulationConfig config = {.ebn0_db=0, .seed=2026, .max_frames=100, .max_frame_errors=100};
    SimulationResult first = {0}, repeat = {0};
    CHECK(simulation_run(&h, &form, &decoder, &config, &first) == DECODE_SUCCESS);
    check_oracle(&h, &config, &first, 0);
    CHECK(first.frames == 100 && first.proven_frames == 100 && !first.truncated_frames);
    CHECK(first.frame_errors > 0 && first.bit_errors >= first.frame_errors);
    CHECK(first.totals.candidates_popped == first.totals.candidates_evaluated);
    CHECK(first.totals.nodes_expanded <= first.totals.candidates_popped && !first.totals.batches);
    CHECK(simulation_run(&h, &form, &decoder, &config, &repeat) == DECODE_SUCCESS);
    char a[2048], b[2048];
    csv_text(&decoder, &config, &first, a, sizeof(a));
    csv_text(&decoder, &config, &repeat, b, sizeof(b));
    CHECK(!strcmp(a, b));
    config.max_frames = 1000;
    config.max_frame_errors = 3;
    CHECK(simulation_run(&h, &form, &decoder, &config, &repeat) == DECODE_SUCCESS);
    check_oracle(&h, &config, &repeat, 0);
    CHECK(repeat.frame_errors == 3 && repeat.frames < 1000 && repeat.error_limit_hit);
    config.max_frames = config.max_frame_errors = 100;
    decoder.candidate_limit = 1;
    CHECK(simulation_run(&h, &form, &decoder, &config, &repeat) == DECODE_SUCCESS);
    check_oracle(&h, &config, &repeat, 1);
    CHECK(repeat.truncated_frames > 0 && repeat.proven_frames + repeat.truncated_frames == 100);
    CHECK(repeat.totals.candidates_evaluated == 100 && repeat.totals.candidates_popped == 100);
    CHECK(repeat.totals.nodes_expanded == 100 && repeat.totals.batches == 0);

    /* Known aggregate: enforce CODEWORD denominator, actual frames and averages. */
    SimulationResult fixture = {.n=8, .k=4, .sigma2=1, .frames=10, .frame_errors=2,
        .bit_errors=4, .hd_bit_errors=8, .proven_frames=10,
        .totals={.candidates_popped=30, .candidates_evaluated=30, .nodes_expanded=15}};
    csv_text(&decoder, &config, &fixture, a, sizeof(a));
    char *row = strchr(a, '\n');
    CHECK(row != NULL);
    char *field = strtok(row + 1, ",\n");
    size_t columns = 0;
    while (field) {
        if (columns == 13) CHECK(strtod(field, NULL) == 0.05);
        if (columns == 14) CHECK(strtod(field, NULL) == 0.2);
        if (columns == 15) CHECK(strtod(field, NULL) == 0.1);
        if (columns == 20) CHECK(strtod(field, NULL) == 3);
        ++columns;
        field = strtok(NULL, ",\n");
    }
    CHECK(columns == 24);
    CHECK(!simulation_write_csv(NULL, &decoder, &config, &first));
    config.max_frames = 0;
    CHECK(simulation_run(&h, &form, &decoder, &config, &repeat) == DECODE_INVALID_ARGUMENT);
    config.max_frames = 100;
    config.max_frame_errors = 0;
    CHECK(simulation_run(&h, &form, &decoder, &config, &repeat) == DECODE_INVALID_ARGUMENT);
    config.max_frame_errors = 100;
    config.ebn0_db = NAN;
    CHECK(simulation_run(&h, &form, &decoder, &config, &repeat) == DECODE_INVALID_ARGUMENT);
    binary_systematic_free(&form);
    printf("Simulation checks passed (%u checks; 100-frame replay, independent ML counts, error cap, truncation and CSV)\n", checks);
    return EXIT_SUCCESS;
}
