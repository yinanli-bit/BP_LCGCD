#include "bplcgcd/channel.h"

#include <math.h>

/* SplitMix64's fixed integer mixing formula; all arithmetic is uint64_t. */
//伪随机数生成
static uint64_t mix64(uint64_t value)
{
    value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31);
}

//生成 (0,1) 之间的均匀随机数
static double uniform_open(uint64_t *state)
{
    *state += UINT64_C(0x9e3779b97f4a7c15);
    /* Midpoints of 2^52 bins stay strictly inside (0,1).
     * Using 53 bits + 0.5 could round the largest value up to 1. */
    return ((double)(mix64(*state) >> 12) + 0.5) * 0x1p-52;
}

ChannelStatus bpsk_noise_variance(double rate, double ebn0_db, double *sigma2)
{
    if (!sigma2 || !isfinite(rate) || rate <= 0 || rate > 1 || !isfinite(ebn0_db)) {
        return CHANNEL_INVALID_ARGUMENT;
    }
    double value = 1.0 / (2.0 * rate * pow(10.0, ebn0_db / 10.0));
    if (!isfinite(value) || value <= 0) return CHANNEL_NUMERIC_ERROR;
    *sigma2 = value;
    return CHANNEL_OK;
}

ChannelStatus bpsk_modulate(const BitVector *codeword, double *symbols)
{
    if (!codeword || !codeword->data || !codeword->length || !symbols) return CHANNEL_INVALID_ARGUMENT;
    for (size_t j = 0; j < codeword->length; ++j) {
        if (codeword->data[j] > 1) return CHANNEL_INVALID_ARGUMENT;
        symbols[j] = 1.0 - 2.0 * codeword->data[j];
    }
    return CHANNEL_OK;
}

ChannelStatus bpsk_awgn_transmit(const BitVector *codeword, double sigma2,
                               uint64_t seed, uint64_t frame_index, double *received)
{
    if (!isfinite(sigma2) || sigma2 < 0) return CHANNEL_INVALID_ARGUMENT;
    ChannelStatus status = bpsk_modulate(codeword, received);
    if (status != CHANNEL_OK || sigma2 == 0) return status;
    uint64_t state = seed ^ mix64(frame_index);
    double sigma = sqrt(sigma2);
    const double two_pi = 6.283185307179586476925286766559;
    for (size_t j = 0; j < codeword->length; ++j) {
        double u = uniform_open(&state), v = uniform_open(&state);
        double gaussian = sqrt(-2.0 * log(u)) * cos(two_pi * v); //Box–Muller 变换
        received[j] += sigma * gaussian;
    }
    return CHANNEL_OK;
}

ChannelStatus bpsk_compute_llr(const double *received, size_t length, double sigma2, double *llr)
{
    if (!received || !llr || !length || !isfinite(sigma2) || sigma2 <= 0) return CHANNEL_INVALID_ARGUMENT;
    for (size_t j = 0; j < length; ++j) {
        if (!isfinite(received[j])) return CHANNEL_INVALID_ARGUMENT;
        double value = 2.0 * (received[j] / sigma2);
        if (!isfinite(value)) return CHANNEL_NUMERIC_ERROR;
        llr[j] = value;
    }
    return CHANNEL_OK;
}
