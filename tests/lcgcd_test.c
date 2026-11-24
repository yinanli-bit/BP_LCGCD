#include "bplcgcd/decoder.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

static unsigned checks, frames;
#define CHECK(condition) do { \
    if (!(condition)) { \
        fprintf(stderr,"%s:%d: %s\n",__FILE__,__LINE__,#condition); \
        exit(EXIT_FAILURE); \
    } \
    ++checks; \
} while (0)

static unsigned syndrome(const BinaryMatrix *h, unsigned bits)
{
    unsigned result=0;
    for (size_t i=0; i<h->rows; ++i) {
        unsigned parity=0;
        for (size_t j=0; j<h->cols; ++j) parity ^= h->data[i*h->cols+j] & ((bits>>j)&1u);
        result |= parity << i;
    }
    return result;
}

static unsigned pack(const BitVector *v)
{
    unsigned result=0;
    for (size_t j=0; j<v->length; ++j) { CHECK(v->data[j]<=1); result |= (unsigned)v->data[j]<<j; }
    return result;
}

/* Full original-coordinate ML oracle, independent of S03/S10/SLVA/completion. */
static void verify(DecodeInput *input, size_t delta, uint64_t limit)
{
    unsigned target=input->mode==DECODE_CLASSICAL ? 0 : pack(input->s_target), reference=0;
    double weights[6], optimum=INFINITY;
    for (size_t j=0; j<input->length; ++j) {
        unsigned bit=input->mode==DECODE_SYNDROME ? (input->z_hint ? input->z_hint->data[j] : 0) : input->llr[j]<0;
        reference |= bit<<j;
        weights[j]=input->mode==DECODE_SYNDROME ? input->weights[j] : fabs(input->llr[j]);
    }
    for (unsigned f=0; f<(1u<<input->length); ++f) if (syndrome(input->h,reference^f)==target) {
        double cost=0;
        for (size_t j=0; j<input->length; ++j) cost += ((f>>j)&1u)*weights[j];
        if (cost<optimum) optimum=cost;
    }
    uint8_t bits[6];
    DecodeResult result={.estimate={input->length,bits}};
    DecodeConfig config={.decoder_name="lcgcd",.delta=delta,.candidate_limit=limit};
    DecodeStatus status=decoder_decode(&config,input,&result);
    CHECK(status==result.status);
    if (!isfinite(optimum)) {
        CHECK(status==DECODE_INCONSISTENT_SYNDROME && !result.candidate_valid && !result.optimality_proven);
        CHECK(!result.stats.candidates_popped);
    } else {
        CHECK(status==DECODE_SUCCESS || status==DECODE_TRUNCATED);
        CHECK(result.candidate_valid && syndrome(input->h,pack(&result.estimate))==target);
        double cost=0;
        for (size_t j=0; j<input->length; ++j) cost += (((pack(&result.estimate)^reference)>>j)&1u)*weights[j];
        CHECK(cost==result.cost && cost>=optimum);
        CHECK(result.stats.candidates_popped==result.stats.candidates_evaluated && result.stats.batches==0);
        CHECK(!limit || result.stats.candidates_evaluated<=limit);
        if (status==DECODE_SUCCESS) CHECK(result.optimality_proven && cost==optimum);
        else CHECK(limit && !result.optimality_proven && result.stop_reason==DECODE_STOP_BUDGET);
        if (!limit) CHECK(status==DECODE_SUCCESS);
    }
    ++frames;
}

static void test_hand_and_modes(void)
{
    uint8_t data[]={0,1,1,0,1,0, 0,0,1,1,0,1, 0,1,0,1,1,1};
    BinaryMatrix h={3,6,data};
    BinarySystematic form={0};
    CHECK(binary_matrix_systematize(&h,&form)==GF2_OK);
    double weights[]={1000,3,6,9,8,7}, llr[]={1000,-3,6,9,-8,7};
    uint8_t target_data[]={1,1,0}, hint_data[]={0,1,0,0,1,0}, bits[6];
    BitVector target={3,target_data}, hint={6,hint_data};
    DecodeInput input={.mode=DECODE_SYNDROME,.h=&h,.length=6,.weights=weights,.s_target=&target,.systematic=&form};
    DecodeResult result={.estimate={6,bits}};
    FILE *trace=tmpfile();
    CHECK(trace);
    DecodeConfig config={.decoder_name="lcgcd",.delta=1,.trace=trace};
    CHECK(decoder_decode(&config,&input,&result)==DECODE_SUCCESS);
    CHECK(result.optimality_proven && result.cost==6 && pack(&result.estimate)==4);
    CHECK(result.stats.candidates_evaluated==1 && result.stop_reason==DECODE_STOP_BOUND);
    rewind(trace);
    char text[4096];
    size_t length=fread(text,1,sizeof(text)-1,trace); text[length]=0;
    CHECK(strstr(text,"f_A/P(stage)=10000 W_local=6 f_I(sys)=01 W=6"));
    CHECK(strstr(text,"trace stop=bound"));
    fclose(trace);
    for (unsigned mode=0; mode<3; ++mode) {
        input.mode=(DecodeMode)mode;
        input.llr=mode==DECODE_SYNDROME ? NULL : llr;
        input.weights=mode==DECODE_SYNDROME ? weights : NULL;
        input.z_hint=mode==DECODE_SYNDROME ? &hint : NULL;
        input.s_target=mode==DECODE_CLASSICAL ? NULL : &target;
        verify(&input,1,0); verify(&input,1,1);
    }
    input=(DecodeInput){.mode=DECODE_SYNDROME,.h=&h,.length=6,.weights=weights,.s_target=&target,.systematic=&form};
    config=(DecodeConfig){.decoder_name="lcgcd",.delta=2,.candidate_limit=1,.slva_max_candidates=1};
    CHECK(decoder_decode(&config,&input,&result)==DECODE_SUCCESS && result.optimality_proven);
    CHECK(result.stats.candidates_evaluated==1 && result.stats.nodes_expanded==0 && result.cost==6);
    config.state_limit=1;
    CHECK(decoder_decode(&config,&input,&result)==DECODE_TRUNCATED && !result.candidate_valid);
    CHECK(result.stop_reason==DECODE_STOP_RESOURCE && !result.optimality_proven);
    config.state_limit=0; config.slva_max_bytes=1;
    CHECK(decoder_decode(&config,&input,&result)==DECODE_TRUNCATED && result.stop_reason==DECODE_STOP_RESOURCE);
    CHECK(!result.candidate_valid);
    config.slva_max_bytes=0; config.delta=1;
    target_data[1]=0; target_data[2]=1; /* Local cost 0, completed Ac cost 3. */
    CHECK(decoder_decode(&config,&input,&result)==DECODE_TRUNCATED && result.candidate_valid);
    CHECK(result.cost==3 && result.stop_reason==DECODE_STOP_RESOURCE && !result.optimality_proven);
    config.slva_max_candidates=0;
    CHECK(decoder_decode(&config,&input,&result)==DECODE_SUCCESS && result.cost==3 && result.optimality_proven);
    config.delta=3;
    CHECK(decoder_decode(&config,&input,&result)==DECODE_DIMENSION_MISMATCH && !result.candidate_valid);
    target_data[2]=0; config=(DecodeConfig){.decoder_name="lcgcd",.delta=1};
    CHECK(decoder_decode(&config,&input,&result)==DECODE_INCONSISTENT_SYNDROME);
    CHECK(!result.candidate_valid && !result.stats.candidates_evaluated);
    binary_systematic_free(&form);
}

static void test_budget(void)
{
    uint8_t data[]={1,0,1,0, 0,1,0,1}, target_data[]={1,1}, bits[4];
    BinaryMatrix h={2,4,data};
    BinarySystematic form={0};
    CHECK(binary_matrix_systematize(&h,&form)==GF2_OK);
    double weights[]={10,10,1,1};
    BitVector target={2,target_data};
    DecodeInput input={.mode=DECODE_SYNDROME,.h=&h,.length=4,.weights=weights,.s_target=&target,.systematic=&form};
    DecodeResult result={.estimate={4,bits}};
    DecodeConfig config={.decoder_name="lcgcd",.delta=1,.candidate_limit=1};
    CHECK(decoder_decode(&config,&input,&result)==DECODE_TRUNCATED && result.cost==11);
    CHECK(result.candidate_valid && !result.optimality_proven && result.stop_reason==DECODE_STOP_BUDGET);
    CHECK(pack(&result.estimate)==6 && result.stats.candidates_evaluated==1);
    config.candidate_limit=2;
    CHECK(decoder_decode(&config,&input,&result)==DECODE_SUCCESS && result.cost==2 && result.optimality_proven);
    for (size_t delta=0; delta<=2; ++delta) { verify(&input,delta,0); verify(&input,delta,1); }
    binary_systematic_free(&form);
}

int main(void)
{
    test_hand_and_modes();
    test_budget();
    puts("LC-GCD hand example: Eq.10 completes Ac; first legal local path alone is not an ML proof");
    printf("LC-GCD checks passed (%u checks, %u focused ML oracle frames)\n",checks,frames);
    return EXIT_SUCCESS;
}
