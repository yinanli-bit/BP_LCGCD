#ifndef BPLCGCD_CHANNEL_H
#define BPLCGCD_CHANNEL_H

#include "bplcgcd/gf2.h"

//码字 c → x=1−2c → y=x+σG → LLR=2y/σ² → ĉ=1[y<0]

typedef enum { 
    CHANNEL_OK = 0, 
    CHANNEL_INVALID_ARGUMENT, 
    CHANNEL_NUMERIC_ERROR 
} ChannelStatus;

//单帧测试/后续仿真驱动：由码率、Eb/N0 得到正方差
ChannelStatus bpsk_noise_variance(double rate, double ebn0_db, double *sigma2);

//Bits/samples stay in original column order; x[j]=1-2*codeword[j].
ChannelStatus bpsk_modulate(const BitVector *codeword, double *symbols);

//y=x+σG
ChannelStatus bpsk_awgn_transmit(const BitVector *codeword, double sigma2,uint64_t seed, uint64_t frame_index, double *received);

// log(P(y|c=0)/P(y|c=1))=2*y/sigma2;
ChannelStatus bpsk_compute_llr(const double *received, size_t length, double sigma2, double *llr);

#endif
