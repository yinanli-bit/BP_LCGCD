#include "bplcgcd/simulation.h"
#include "bplcgcd/channel.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

static unsigned checks, frames_checked;
#define CHECK(condition) do { if (!(condition)) { \
    fprintf(stderr,"%s:%d: %s\n",__FILE__,__LINE__,#condition); exit(EXIT_FAILURE); \
} ++checks; } while (0)

static void same_summary(const SimulationResult *a, const SimulationResult *b)
{
    CHECK(a->n==b->n && a->k==b->k && a->sigma2==b->sigma2 && a->frames==b->frames);
    CHECK(a->frame_errors==b->frame_errors && a->bit_errors==b->bit_errors && a->hd_bit_errors==b->hd_bit_errors);
    CHECK(a->proven_frames==b->proven_frames && a->truncated_frames==b->truncated_frames);
    CHECK(a->error_limit_hit==b->error_limit_hit);
    CHECK(a->totals.candidates_popped==b->totals.candidates_popped);
    CHECK(a->totals.candidates_evaluated==b->totals.candidates_evaluated);
    CHECK(a->totals.nodes_expanded==b->totals.nodes_expanded && a->totals.batches==b->totals.batches);
    CHECK(a->totals.partial_patterns_generated==b->totals.partial_patterns_generated);
    CHECK(a->totals.trellis_states_visited==b->totals.trellis_states_visited);
    CHECK(a->p50_queries==b->p50_queries && a->p95_queries==b->p95_queries && a->p99_queries==b->p99_queries);
    CHECK(!memcmp(a->early_stop_counts,b->early_stop_counts,sizeof(a->early_stop_counts)));
}

/* Independent empirical-CDF definition, without the production sort/rank code. */
static uint64_t quantile(const uint64_t *values, size_t count, unsigned percent)
{
    uint64_t best=UINT64_MAX;
    for (size_t i=0; i<count; ++i) {
        size_t below=0;
        for (size_t j=0; j<count; ++j)
            below+=values[j]<=values[i];
        if (100*below>=percent*count && values[i]<best)
            best=values[i];
    }
    return best;
}

static unsigned pack(const BitVector *v)
{
    unsigned word=0;
    for (size_t j=0; j<v->length; ++j) { CHECK(v->data[j]<=1); word |= (unsigned)v->data[j]<<j; }
    return word;
}

/* Independent original-matrix parity and Euclidean-distance ML oracle. */
static unsigned ml_word(const BinaryMatrix *h, const double *y)
{
    unsigned chosen=0;
    double best=INFINITY;
    for (unsigned word=0; word<128; ++word) {
        unsigned invalid=0;
        for (size_t i=0; i<h->rows; ++i) {
            unsigned parity=0;
            for (size_t j=0; j<7; ++j) parity ^= h->data[i*7+j]&((word>>j)&1u);
            invalid |= parity;
        }
        if (invalid) continue;
        double distance=0;
        for (size_t j=0; j<7; ++j) { double e=y[j]-(1-2*(int)((word>>j)&1u)); distance+=e*e; }
        if (distance<best) { best=distance; chosen=word; }
    }
    return chosen;
}

static void samples(const char *line, const char *label, const double *expected)
{
    const char *p=strstr(line,label);
    CHECK(p); p+=strlen(label);
    for (size_t j=0; j<7; ++j) {
        char *end;
        double value=strtod(p,&end);
        CHECK(end!=p && !memcmp(&value,&expected[j],sizeof(value)));
        p=end;
        if (j<6) { CHECK(*p==','); ++p; }
    }
}

static void same_text(FILE *a, FILE *b)
{
    rewind(a); rewind(b);
    char x[4096], y[4096];
    while (fgets(x,sizeof(x),a)) { CHECK(fgets(y,sizeof(y),b)); CHECK(!strcmp(x,y)); }
    CHECK(feof(a) && fgetc(b)==EOF);
}

static void check_frames(const BinaryMatrix *h, const BinarySystematic *form,
                         const DecodeConfig *lc, const SimulationConfig *config,
                         const SimulationResult actual[2], FILE *trace)
{
    CHECK(actual[0].frames==actual[1].frames && actual[0].frames>0);
    double sigma2=1/(2*(4.0/7)*pow(10,config->ebn0_db/10));
    SimulationResult expected[2]={{.n=7,.k=4,.sigma2=sigma2},{.n=7,.k=4,.sigma2=sigma2}};
    uint8_t zero_data[7]={0}, bits[2][7];
    BitVector zero={7,zero_data};
    DecodeResult decoded[2]={{.estimate={7,bits[0]}},{.estimate={7,bits[1]}}};
    DecodeConfig decoders[2]={*lc,*lc}; decoders[0].decoder_name="gcd";
    double y[7], llr[7], noise[7];
    uint64_t queries[2][100];
    char input_lines[2][4096], result_lines[2][4096], line[4096];
    rewind(trace);
    CHECK(fgets(line,sizeof(line),trace) && strstr(line,"paired_config "));
    for (uint64_t frame=0; frame<actual[0].frames; ++frame) {
        CHECK(expected[0].frame_errors<config->max_frame_errors && expected[1].frame_errors<config->max_frame_errors);
        CHECK(bpsk_awgn_transmit(&zero,sigma2,config->seed,frame,y)==CHANNEL_OK);
        CHECK(bpsk_compute_llr(y,7,sigma2,llr)==CHANNEL_OK);
        for (size_t j=0; j<7; ++j) noise[j]=y[j]-1;
        DecodeInput in={.mode=DECODE_CLASSICAL,.h=h,.length=7,.llr=llr,.systematic=form};
        for (size_t d=0; d<2; ++d) {
            CHECK(fgets(input_lines[d],sizeof(input_lines[d]),trace));
            char prefix[100];
            snprintf(prefix,sizeof(prefix),"paired_input frame=%llu decoder=%s tx=0000000 ",
                     (unsigned long long)frame,decoders[d].decoder_name);
            CHECK(!strncmp(input_lines[d],prefix,strlen(prefix)));
            DecodeStatus status=decoder_decode(&decoders[d],&in,&decoded[d]);
            CHECK((status==DECODE_SUCCESS || status==DECODE_TRUNCATED) && decoded[d].candidate_valid);
            CHECK(!lc->candidate_limit || decoded[d].stats.candidates_evaluated<=lc->candidate_limit);
        }
        CHECK(!strcmp(strstr(input_lines[0]," tx="),strstr(input_lines[1]," tx=")));
        samples(input_lines[0]," noise=",noise); samples(input_lines[0]," y=",y); samples(input_lines[0]," llr=",llr);
        unsigned optimum=ml_word(h,y);
        for (size_t d=0; d<2; ++d) {
            CHECK(fgets(result_lines[d],sizeof(result_lines[d]),trace));
            char prefix[100], estimate[8];
            snprintf(prefix,sizeof(prefix),"paired_result frame=%llu decoder=%s status=%s estimate=",
                     (unsigned long long)frame,decoders[d].decoder_name,decoded[d].status==DECODE_SUCCESS ? "SUCCESS" : "TRUNCATED");
            CHECK(!strncmp(result_lines[d],prefix,strlen(prefix)));
            for (size_t j=0; j<7; ++j) estimate[j]=(char)('0'+decoded[d].estimate.data[j]);
            estimate[7]=0;
            CHECK(!strncmp(result_lines[d]+strlen(prefix),estimate,7));
            char *cost=strstr(result_lines[d]," cost="); CHECK(cost);
            CHECK(strtod(cost+6,NULL)==decoded[d].cost);
            CHECK(fgets(line,sizeof(line),trace) && strstr(line,"paired_stats "));
            unsigned word=pack(&decoded[d].estimate), errors=0;
            if (decoded[d].optimality_proven) CHECK(word==optimum);
            for (unsigned copy=word; copy; copy>>=1) errors+=copy&1u;
            ++expected[d].frames; expected[d].bit_errors+=errors; expected[d].frame_errors+=word!=0;
            for (size_t j=0; j<7; ++j) expected[d].hd_bit_errors+=y[j]<0;
            expected[d].proven_frames+=decoded[d].optimality_proven;
            expected[d].truncated_frames+=decoded[d].status==DECODE_TRUNCATED;
            expected[d].totals.candidates_popped+=decoded[d].stats.candidates_popped;
            expected[d].totals.candidates_evaluated+=decoded[d].stats.candidates_evaluated;
            expected[d].totals.nodes_expanded+=decoded[d].stats.nodes_expanded;
            expected[d].totals.batches+=decoded[d].stats.batches;
            expected[d].totals.partial_patterns_generated+=decoded[d].stats.partial_patterns_generated;
            expected[d].totals.trellis_states_visited+=decoded[d].stats.trellis_states_visited;
            ++expected[d].early_stop_counts[decoded[d].stop_reason];
            queries[d][frame]=decoded[d].stats.candidates_evaluated;
        }
        if (!lc->delta) CHECK(!strcmp(strstr(result_lines[0]," status="),strstr(result_lines[1]," status=")));
        ++frames_checked;
    }
    CHECK(fgetc(trace)==EOF);
    bool stop=expected[0].frame_errors>=config->max_frame_errors || expected[1].frame_errors>=config->max_frame_errors;
    CHECK(actual[0].frames==config->max_frames || stop);
    for (size_t d=0; d<2; ++d) {
        expected[d].error_limit_hit=stop;
        expected[d].p50_queries=quantile(queries[d],actual[d].frames,50);
        expected[d].p95_queries=quantile(queries[d],actual[d].frames,95);
        expected[d].p99_queries=quantile(queries[d],actual[d].frames,99);
        same_summary(&actual[d],&expected[d]);
        CHECK(actual[d].runtime_available && isfinite(actual[d].runtime_seconds) && actual[d].runtime_seconds>=0);
        CHECK(actual[d].proven_frames==actual[d].early_stop_counts[DECODE_STOP_BOUND]+actual[d].early_stop_counts[DECODE_STOP_EXHAUSTED]);
        CHECK(actual[d].truncated_frames==actual[d].early_stop_counts[DECODE_STOP_BUDGET]+actual[d].early_stop_counts[DECODE_STOP_RESOURCE]);
    }
}

static void check_csv(const DecodeConfig *lc, const SimulationConfig *config, const SimulationResult result[2])
{
    FILE *csv=tmpfile(); CHECK(csv);
    CHECK(simulation_compare_write_csv(csv,lc,config,result));
    rewind(csv);
    char header[2048], rows[2][2048], prefix[40];
    CHECK(fgets(header,sizeof(header),csv) && !strncmp(header,"schema_version,decoder,delta,",29));
    for (size_t d=0; d<2; ++d) {
        CHECK(fgets(rows[d],sizeof(rows[d]),csv));
        snprintf(prefix,sizeof(prefix),"1,%s,%zu,",d ? "lcgcd" : "gcd",d ? lc->delta : 0);
        CHECK(!strncmp(rows[d],prefix,strlen(prefix)));
        unsigned columns=1;
        for (char *p=rows[d]; *p; ++p) columns+=*p==',';
        CHECK(columns==40);
        const char *runtime=strstr(header,"runtime_seconds");
        CHECK(runtime);
        size_t column=0;
        for (const char *p=header; p<runtime; ++p)
            column+=*p==',';
        char *value=rows[d];
        for (size_t i=0; i<column; ++i)
            value=strchr(value,',')+1;
        char *end=strchr(value,',');
        CHECK(end);
        memmove(value,end,strlen(end)+1);
    }
    CHECK(fgetc(csv)==EOF);
    if (!lc->delta) {
        char *a=strchr(strchr(rows[0],',')+1,',');
        char *b=strchr(strchr(rows[1],',')+1,',');
        CHECK(!strcmp(a,b));
    }
    fclose(csv);
}

int main(int argc, char **argv)
{
    if (argc>1) CHECK(freopen(argv[1],"w",stdout));
    uint8_t entries[]={1,0,0,1,1,0,1, 0,1,0,1,0,1,1,
                       0,0,1,0,1,1,1, 1,1,0,0,1,1,0};
    BinaryMatrix h={4,7,entries};
    BinarySystematic form={0};
    CHECK(binary_matrix_systematize(&h,&form)==GF2_OK && form.rank==3);
    SimulationConfig config={.seed=2026,.ebn0_db=0,.max_frames=12,.max_frame_errors=12};
    DecodeConfig lc={.decoder_name="lcgcd"};
    SimulationResult result[2], repeat[2], single;
    for (size_t delta=0; delta<=3; ++delta) for (unsigned budget=0; budget<2; ++budget) {
        lc.delta=delta; lc.candidate_limit=budget;
        FILE *trace=tmpfile(), *again=tmpfile(); CHECK(trace && again);
        CHECK(simulation_compare(&h,&form,&lc,&config,result,trace)==DECODE_SUCCESS);
        check_frames(&h,&form,&lc,&config,result,trace); check_csv(&lc,&config,result);
        CHECK(simulation_compare(&h,&form,&lc,&config,repeat,again)==DECODE_SUCCESS);
        for (size_t d=0; d<2; ++d) same_summary(&result[d],&repeat[d]);
        same_text(trace,again); fclose(trace); fclose(again);
        if (!delta) {
            same_summary(&result[0],&result[1]);
            DecodeConfig gcd=lc; gcd.decoder_name="gcd";
            CHECK(simulation_run(&h,&form,&gcd,&config,&single)==DECODE_SUCCESS);
            same_summary(&result[0],&single);
        }
    }
    lc.delta=1;
    lc.candidate_limit=0;
    lc.slva_max_candidates=1;
    FILE *resource_trace=tmpfile();
    CHECK(resource_trace);
    CHECK(simulation_compare(&h,&form,&lc,&config,result,resource_trace)==DECODE_SUCCESS);
    check_frames(&h,&form,&lc,&config,result,resource_trace);
    CHECK(result[1].early_stop_counts[DECODE_STOP_RESOURCE]>0 && !result[1].early_stop_counts[DECODE_STOP_BUDGET]);
    fclose(resource_trace);
    lc.slva_max_candidates=0;
    lc.candidate_limit=1; config.ebn0_db=-3; config.max_frames=100; config.max_frame_errors=2;
    bool asymmetric=false;
    for (unsigned trial=0; trial<8 && !asymmetric; ++trial) {
        config.seed=2026+trial;
        FILE *trace=tmpfile(); CHECK(trace);
        CHECK(simulation_compare(&h,&form,&lc,&config,result,trace)==DECODE_SUCCESS);
        check_frames(&h,&form,&lc,&config,result,trace); fclose(trace);
        CHECK(result[0].error_limit_hit && result[1].error_limit_hit && result[0].frames<100);
        asymmetric=result[0].frame_errors!=result[1].frame_errors;
    }
    CHECK(asymmetric);
    lc.delta=2; lc.state_limit=1;
    CHECK(simulation_compare(&h,&form,&lc,&config,result,NULL)==DECODE_INVALID_CANDIDATE);
    CHECK(!result[0].frames && !result[1].frames && !result[0].totals.candidates_popped);
    CHECK(!simulation_compare_write_csv(NULL,&lc,&config,result));
    CHECK(simulation_compare(&h,&form,NULL,&config,result,NULL)==DECODE_INVALID_ARGUMENT);
    binary_systematic_free(&form);
    fprintf(stderr,"Paired simulation passed (%u checks, %u independently replayed frame pairs; shared samples, delta=0, caps and CSV)\n",checks,frames_checked);
    return EXIT_SUCCESS;
}
