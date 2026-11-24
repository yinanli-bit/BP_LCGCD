#include "bplcgcd/channel.h"
#include "bplcgcd/decoder.h"

#include <float.h>
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

static int near(double a, double b)
{
    return fabs(a - b) <= 1e-12 * fmax(1.0, fabs(b));
}

static void test_parameters_and_llr(void)
{
    double variance = -1;
    CHECK(bpsk_noise_variance(0.5, 0, &variance) == CHANNEL_OK && variance == 1);
    CHECK(bpsk_noise_variance(1, 0, &variance) == CHANNEL_OK && variance == 0.5);
    CHECK(bpsk_noise_variance(4.0 / 7, 0, &variance) == CHANNEL_OK && variance == 0.875);
    CHECK(bpsk_noise_variance(0.5, 10, &variance) == CHANNEL_OK && near(variance, 0.1));
    CHECK(bpsk_noise_variance(0.5, -10, &variance) == CHANNEL_OK && near(variance, 10));
    CHECK(bpsk_noise_variance(0.5, 0, NULL) == CHANNEL_INVALID_ARGUMENT);
    const double bad_rates[] = {0, -1, 1.1, NAN, INFINITY};
    for (size_t i = 0; i < sizeof(bad_rates) / sizeof(*bad_rates); ++i) {
        CHECK(bpsk_noise_variance(bad_rates[i], 0, &variance) == CHANNEL_INVALID_ARGUMENT);
    }
    CHECK(bpsk_noise_variance(0.5, NAN, &variance) == CHANNEL_INVALID_ARGUMENT);
    CHECK(bpsk_noise_variance(0.5, INFINITY, &variance) == CHANNEL_INVALID_ARGUMENT);
    CHECK(bpsk_noise_variance(0.5, 4000, &variance) == CHANNEL_NUMERIC_ERROR);
    CHECK(bpsk_noise_variance(0.5, -4000, &variance) == CHANNEL_NUMERIC_ERROR);

    uint8_t data[] = {0, 1, 0, 1};
    BitVector word = {4, data};
    double x[4], y[4], llr[4];
    CHECK(bpsk_modulate(&word, x) == CHANNEL_OK);
    CHECK(x[0] == 1 && x[1] == -1 && x[2] == 1 && x[3] == -1);
    CHECK(bpsk_awgn_transmit(&word, 0, 2026, 7, y) == CHANNEL_OK);
    CHECK(memcmp(x, y, sizeof(x)) == 0);
    CHECK(bpsk_compute_llr(y, 4, 0.5, llr) == CHANNEL_OK);
    CHECK(llr[0] == 4 && llr[1] == -4 && llr[2] == 4 && llr[3] == -4);
    /* Independent likelihood-ratio definition, including y=0 hard-decision tie. */
    const double received[] = {1.25, -0.75, 0, -0.0};
    CHECK(bpsk_compute_llr(received, 4, 0.875, llr) == CHANNEL_OK);
    for (size_t j = 0; j < 4; ++j) {
        double yj = received[j];
        double p0 = exp(-(yj - 1) * (yj - 1) / 1.75);
        double p1 = exp(-(yj + 1) * (yj + 1) / 1.75);
        CHECK(near(llr[j], log(p0 / p1)));
        CHECK((llr[j] < 0) == (yj < 0));
    }
    CHECK(bpsk_modulate(NULL, x) == CHANNEL_INVALID_ARGUMENT);
    CHECK(bpsk_modulate(&word, NULL) == CHANNEL_INVALID_ARGUMENT);
    BitVector empty = {0};
    CHECK(bpsk_modulate(&empty, x) == CHANNEL_INVALID_ARGUMENT);
    data[3] = 2;
    CHECK(bpsk_modulate(&word, x) == CHANNEL_INVALID_ARGUMENT);
    data[3] = 1;
    const double invalid_variances[] = {-1, NAN, INFINITY};
    for (size_t i = 0; i < 3; ++i) {
        CHECK(bpsk_awgn_transmit(&word, invalid_variances[i], 0, 0, y) == CHANNEL_INVALID_ARGUMENT);
        CHECK(bpsk_compute_llr(y, 4, invalid_variances[i], llr) == CHANNEL_INVALID_ARGUMENT);
    }
    CHECK(bpsk_compute_llr(y, 4, 0, llr) == CHANNEL_INVALID_ARGUMENT);
    CHECK(bpsk_compute_llr(NULL, 4, 1, llr) == CHANNEL_INVALID_ARGUMENT);
    CHECK(bpsk_compute_llr(y, 4, 1, NULL) == CHANNEL_INVALID_ARGUMENT);
    CHECK(bpsk_compute_llr(y, 0, 1, llr) == CHANNEL_INVALID_ARGUMENT);
    const double nonfinite = NAN, huge = DBL_MAX;
    CHECK(bpsk_compute_llr(&nonfinite, 1, 1, llr) == CHANNEL_INVALID_ARGUMENT);
    CHECK(bpsk_compute_llr(&huge, 1, 1, llr) == CHANNEL_NUMERIC_ERROR);
    CHECK(bpsk_compute_llr(&huge, 1, DBL_MAX, llr) == CHANNEL_OK && llr[0] == 2);
}

static void test_replay_and_noise(void)
{
    uint8_t bits[7] = {0};
    BitVector word = {7, bits};
    double a[7], b[7], other[7];
    CHECK(bpsk_awgn_transmit(&word, 1, 0, 0, a) == CHANNEL_OK);
    /* Recorded seed=0/frame=0 samples freeze RNG order and Gaussian mapping. */
    const double first[] = {0.54724225978254171, 3.6506058120796689,
                            0.011395875375671549, 1.2524627241506141};
    for (size_t j = 0; j < 4; ++j) CHECK(near(a[j], first[j]));
    CHECK(bpsk_awgn_transmit(&word, 1, 0, 1, other) == CHANNEL_OK);
    CHECK(memcmp(a, other, sizeof(a)) != 0);
    CHECK(bpsk_awgn_transmit(&word, 1, 0, 0, b) == CHANNEL_OK);
    CHECK(memcmp(a, b, sizeof(a)) == 0); /* Independent of intervening frame calls. */
    CHECK(bpsk_awgn_transmit(&word, 1, 1, 0, other) == CHANNEL_OK);
    CHECK(memcmp(a, other, sizeof(a)) != 0);
    CHECK(bpsk_awgn_transmit(&word, 4, 0, 0, b) == CHANNEL_OK);
    for (size_t j = 0; j < 7; ++j) CHECK(near(b[j] - 1, 2 * (a[j] - 1)));
    memset(bits, 1, sizeof(bits));
    CHECK(bpsk_awgn_transmit(&word, 1, 0, 0, b) == CHANNEL_OK);
    for (size_t j = 0; j < 7; ++j) CHECK(near(b[j] + 1, a[j] - 1));
    for (size_t j = 0; j < 7; ++j) CHECK(bits[j] == 1);
    CHECK(bpsk_awgn_transmit(&word, 1, UINT64_MAX, UINT64_MAX, a) == CHANNEL_OK);
    CHECK(bpsk_awgn_transmit(&word, 1, UINT64_MAX, UINT64_MAX, b) == CHANNEL_OK);
    CHECK(memcmp(a, b, sizeof(a)) == 0);

    /* One fixed sampler frame; a quick distribution sanity check, not FER. */
    const size_t length = 8192;
    uint8_t *zeros = calloc(length, sizeof(*zeros));
    double *received = malloc(length * sizeof(*received));
    CHECK(zeros && received);
    BitVector zero_word = {length, zeros};
    CHECK(bpsk_awgn_transmit(&zero_word, 0.875, 2026, 23, received) == CHANNEL_OK);
    double sum = 0, squares = 0;
    for (size_t j = 0; j < length; ++j) {
        double g = (received[j] - 1) / sqrt(0.875);
        sum += g;
        squares += g * g;
    }
    double mean = sum / length, variance = squares / length - mean * mean;
    CHECK(fabs(mean) < 0.06 && fabs(variance - 1) < 0.08);
    printf("Fixed noise frame: normalized mean=%.6f variance=%.6f (8192 samples)\n", mean, variance);
    free(received);
    free(zeros);
}

static double distance(const double *received, const uint8_t *word, size_t length)
{
    double value = 0;
    for (size_t j = 0; j < length; ++j) {
        double delta = received[j] - (1.0 - 2.0 * word[j]);
        value += delta * delta;
    }
    return value;
}

static void test_single_frame_gcd(void)
{
    uint8_t entries[] = {1,0,0,1,1,0,1, 0,1,0,1,0,1,1, 0,0,1,0,1,1,1};
    BinaryMatrix h = {3, 7, entries};
    BinarySystematic form = {0};
    CHECK(binary_matrix_systematize(&h, &form) == GF2_OK);
    double sigma2, received[7], llr[7];
    CHECK(bpsk_noise_variance((7.0 - form.rank) / 7, 0, &sigma2) == CHANNEL_OK);
    CHECK(sigma2 == 0.875);
    uint8_t words[2][7] = {{0}, {1,0,1,0,1,0,0}}, decoded[7];
    DecodeInput input = {.mode=DECODE_CLASSICAL, .h=&h, .length=7, .llr=llr, .systematic=&form};
    DecodeConfig config = {.decoder_name="gcd"};
    DecodeResult result = {.estimate={7, decoded}};
    for (size_t i = 0; i < 2; ++i) {
        BitVector word = {7, words[i]};
        CHECK(bpsk_awgn_transmit(&word, 0, 2026, i, received) == CHANNEL_OK);
        CHECK(bpsk_compute_llr(received, 7, sigma2, llr) == CHANNEL_OK);
        CHECK(decoder_decode(&config, &input, &result) == DECODE_SUCCESS);
        CHECK(result.candidate_valid && result.optimality_proven && result.cost == 0);
        CHECK(memcmp(decoded, word.data, 7) == 0);
    }
    BitVector word = {7, words[1]};
    CHECK(bpsk_awgn_transmit(&word, sigma2, 2026, 8, received) == CHANNEL_OK);
    CHECK(bpsk_compute_llr(received, 7, sigma2, llr) == CHANNEL_OK);
    CHECK(decoder_decode(&config, &input, &result) == DECODE_SUCCESS);
    CHECK(result.candidate_valid && result.optimality_proven);
    /* Independent AWGN ML oracle: squared Euclidean distance over codewords,
     * without production LLR, flip scoring, systematization or TEP ordering. */
    double optimum = INFINITY;
    for (unsigned bits = 0; bits < 128; ++bits) {
        uint8_t candidate[7];
        for (size_t j = 0; j < 7; ++j) candidate[j] = (uint8_t)((bits >> j) & 1u);
        unsigned invalid = 0;
        for (size_t i = 0; i < 3; ++i) {
            unsigned parity = 0;
            for (size_t j = 0; j < 7; ++j) parity ^= h.data[i * 7 + j] & candidate[j];
            invalid |= parity;
        }
        if (!invalid) {
            double metric = distance(received, candidate, 7);
            if (metric < optimum) optimum = metric;
        }
    }
    CHECK(near(distance(received, decoded, 7), optimum));
    printf("Single-frame Eb/N0=0 dB, R=4/7: sigma2=%.3f, y[0]=%.6f, LLR[0]=%.6f; ML oracle passed\n",
           sigma2, received[0], llr[0]);
    binary_systematic_free(&form);
}

int main(void)
{
    test_parameters_and_llr();
    test_replay_and_noise();
    test_single_frame_gcd();
    printf("Channel checks passed (%u checks; noiseless, replay, noise scale and AWGN ML oracle)\n", checks);
    return EXIT_SUCCESS;
}
