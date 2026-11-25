#include "bplcgcd/simulation.h"
#include "bplcgcd/channel.h"

#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static const char *const stop_names[]={"none","bound","exhausted","budget","resource"};

static int compare_queries(const void *left, const void *right)
{
    uint64_t a=*(const uint64_t *)left, b=*(const uint64_t *)right;
    return (a>b)-(a<b);
}

static uint64_t percentile(const uint64_t *sorted, size_t count, unsigned percent)
{
    size_t rank=(count/100)*percent+((count%100)*percent+99)/100;
    return sorted[rank-1];
}

static void trace_input(FILE *stream, uint64_t frame, const char *decoder,
                        const BitVector *sent, const double *y, const double *llr)
{
    fprintf(stream,"paired_input frame=%" PRIu64 " decoder=%s tx=",frame,decoder);
    for (size_t j=0; j<sent->length; ++j) 
        fputc('0'+sent->data[j],stream);
    fputs(" noise=",stream);
    for (size_t j=0; j<sent->length; ++j) 
        fprintf(stream,"%s%a",j ? "," : "",y[j]-(1-2*sent->data[j]));
    fputs(" y=",stream);
    for (size_t j=0; j<sent->length; ++j) 
        fprintf(stream,"%s%a",j ? "," : "",y[j]);
    fputs(" llr=",stream);
    for (size_t j=0; j<sent->length; ++j) 
        fprintf(stream,"%s%a",j ? "," : "",llr[j]);
    fputc('\n',stream);
}

static void trace_result(FILE *stream, uint64_t frame, const char *decoder, const DecodeResult *result)
{
    fprintf(stream,"paired_result frame=%" PRIu64 " decoder=%s status=%s estimate=",frame,decoder,
            result->status==DECODE_SUCCESS ? "SUCCESS" : "TRUNCATED");
    for (size_t j=0; j<result->estimate.length; ++j) fputc('0'+result->estimate.data[j],stream);
    fprintf(stream," cost=%a valid=1 proven=%d stop=%s popped=%" PRIu64 " evaluated=%" PRIu64
            " expanded=%" PRIu64 " batches=%" PRIu64 "\n",result->cost,(int)result->optimality_proven,
            stop_names[result->stop_reason],result->stats.candidates_popped,result->stats.candidates_evaluated,
            result->stats.nodes_expanded,result->stats.batches);
    fprintf(stream,"paired_stats frame=%" PRIu64 " decoder=%s partial_patterns_generated=%" PRIu64
            " full_candidates_evaluated=%" PRIu64 " query_batches=%" PRIu64 " trellis_states_visited=%" PRIu64
            " early_stop_reason=%s query_limit_hit=%d\n",frame,decoder,result->stats.partial_patterns_generated,
            result->stats.candidates_evaluated,result->stats.batches,result->stats.trellis_states_visited,
            stop_names[result->stop_reason],result->stop_reason==DECODE_STOP_BUDGET);
}

/* count is 1 or 2, supplied only by the two public entry points below. */
static DecodeStatus run_frames(const BinaryMatrix *h, const BinarySystematic *form,
                               const DecodeConfig *decoders, size_t count, const SimulationConfig *config,
                               SimulationResult *results, FILE *paired_trace)
{
    if (!h || !h->data || !h->cols || !form || form->rank >= h->cols || !decoders ||!config || !results || !config->max_frames ||!config->max_frame_errors || h->cols > SIZE_MAX / sizeof(double) ||config->max_frames > UINT64_MAX / h->cols) 
        return DECODE_INVALID_ARGUMENT;
    if (form->reduced.rows != h->rows || form->reduced.cols != h->cols) 
        return DECODE_DIMENSION_MISMATCH;
    if (config->max_frames > SIZE_MAX/sizeof(uint64_t))
        return DECODE_INVALID_ARGUMENT;
    size_t n = h->cols, k = n - form->rank;
    double sigma2;
    ChannelStatus channel = bpsk_noise_variance((double)k / n, config->ebn0_db, &sigma2);
    if (channel != CHANNEL_OK) 
        return channel == CHANNEL_NUMERIC_ERROR ? DECODE_NUMERIC_ERROR : DECODE_INVALID_ARGUMENT;
    for (size_t d=0; d<count; ++d) 
        results[d]=(SimulationResult){.n=n,.k=k,.sigma2=sigma2,.runtime_available=true};
    BitVector zero = {0};
    DecodeResult decoded[2] = {0};
    double *received = malloc(n * sizeof(*received)), *llr = malloc(n * sizeof(*llr));
    uint64_t *queries[2]={NULL,NULL};
    DecodeStatus status = DECODE_ALLOCATION_FAILED;
    if (!received || !llr || bit_vector_init(&zero, n) != GF2_OK) 
        goto finish;
    for (size_t d=0; d<count; ++d) {
        queries[d]=malloc((size_t)config->max_frames*sizeof(*queries[d]));
        if (!queries[d] || bit_vector_init(&decoded[d].estimate,n)!=GF2_OK)
            goto finish;
    }
    DecodeInput input = {.mode=DECODE_CLASSICAL, .h=h, .length=n, .llr=llr, .systematic=form};
    if (paired_trace) 
        fprintf(paired_trace,"paired_config seed=%" PRIu64 " ebn0_db=%.17g n=%zu k=%zu rate=%a sigma2=%a delta=%zu candidate_limit=%" PRIu64 "\n",
        config->seed,config->ebn0_db,n,k,(double)k/n,sigma2,decoders[1].delta,decoders[1].candidate_limit);
    while (results[0].frames < config->max_frames) {
        uint64_t frame=results[0].frames;
        double elapsed[2]={0,0};
        bool timed[2]={false,false};
        channel = bpsk_awgn_transmit(&zero, sigma2, config->seed, frame, received);
        if (channel == CHANNEL_OK) 
            channel = bpsk_compute_llr(received, n, sigma2, llr);
        if (channel != CHANNEL_OK) 
            { status = DECODE_NUMERIC_ERROR; goto finish; }
        for (size_t d=0; d<count; ++d) {
            if (paired_trace) 
                trace_input(paired_trace,frame,decoders[d].decoder_name,&zero,received,llr);
            clock_t begin=clock();
            status=decoder_decode(&decoders[d],&input,&decoded[d]);
            clock_t end=clock();
            timed[d]=begin!=(clock_t)-1 && end!=(clock_t)-1 && end>=begin;
            if (timed[d])
                elapsed[d]=((double)end-(double)begin)/CLOCKS_PER_SEC;
            if (status!=DECODE_SUCCESS && status!=DECODE_TRUNCATED) 
                goto finish;
            if (!decoded[d].candidate_valid) 
                { status=DECODE_INVALID_CANDIDATE; goto finish; }
        }
        bool stop=false;
        /* Commit the common frame only after BOTH decoders succeeded. */
        for (size_t d=0; d<count; ++d) {
            SimulationResult *result=&results[d];
            queries[d][result->frames]=decoded[d].stats.candidates_evaluated;
            uint64_t errors=0;
            for (size_t j=0; j<n; ++j) {
                errors+=decoded[d].estimate.data[j]; /* Compare against transmitted zero. */
                result->hd_bit_errors+=llr[j]<0;
            }
            ++result->frames;
            result->bit_errors+=errors;
            result->frame_errors+=errors!=0;
            result->truncated_frames+=decoded[d].status==DECODE_TRUNCATED;
            result->proven_frames+=decoded[d].optimality_proven;
            result->totals.candidates_popped+=decoded[d].stats.candidates_popped;
            result->totals.candidates_evaluated+=decoded[d].stats.candidates_evaluated;
            result->totals.nodes_expanded+=decoded[d].stats.nodes_expanded;
            result->totals.batches+=decoded[d].stats.batches;
            result->totals.partial_patterns_generated+=decoded[d].stats.partial_patterns_generated;
            result->totals.trellis_states_visited+=decoded[d].stats.trellis_states_visited;
            ++result->early_stop_counts[decoded[d].stop_reason];
            result->runtime_available &= timed[d];
            result->runtime_seconds+=elapsed[d];
            stop |= result->frame_errors>=config->max_frame_errors;
            if (paired_trace) trace_result(paired_trace,frame,decoders[d].decoder_name,&decoded[d]);
        }
        if (stop) {
            for (size_t d=0; d<count; ++d) 
                results[d].error_limit_hit=true;
            break;
        }
    }
    for (size_t d=0; d<count; ++d) {
        size_t frames=(size_t)results[d].frames;
        qsort(queries[d],frames,sizeof(*queries[d]),compare_queries);
        results[d].p50_queries=percentile(queries[d],frames,50);
        results[d].p95_queries=percentile(queries[d],frames,95);
        results[d].p99_queries=percentile(queries[d],frames,99);
    }
    status = DECODE_SUCCESS; /* The run completed; proof/truncation remain per frame. */
finish:
    for (size_t d=0; d<count; ++d) {
        bit_vector_free(&decoded[d].estimate);
        free(queries[d]);
    }
    bit_vector_free(&zero);
    free(llr);
    free(received);
    return status;
}

DecodeStatus simulation_run(const BinaryMatrix *h, const BinarySystematic *form,
                            const DecodeConfig *decoder, const SimulationConfig *config, SimulationResult *result)
{
    if (!decoder || !decoder->decoder_name) 
        return DECODE_INVALID_ARGUMENT;
    if (strcmp(decoder->decoder_name,"gcd")) 
        return DECODE_NOT_IMPLEMENTED;
    return run_frames(h,form,decoder,1,config,result,NULL);
}

DecodeStatus simulation_compare(const BinaryMatrix *h, const BinarySystematic *form,
                               const DecodeConfig *lc_config, const SimulationConfig *config,
                               SimulationResult results[2], FILE *paired_trace)
{
    if (!lc_config || !lc_config->decoder_name || strcmp(lc_config->decoder_name,"lcgcd")) 
        return DECODE_INVALID_ARGUMENT;
    DecodeConfig decoders[2]={*lc_config,*lc_config};
    decoders[0].decoder_name="gcd";
    return run_frames(h,form,decoders,2,config,results,paired_trace);
}

static const char csv_header[]="schema_version,decoder,delta,ebn0_db,seed,n,k,rate,sigma2,"
    "max_frames,max_frame_errors,candidate_limit,frames,frame_errors,bit_errors,hd_bit_errors,fer,ber,hd_ber,"
    "partial_patterns_generated,full_candidates_evaluated,query_batches,trellis_states_visited,"
    "patterns_popped,nodes_expanded,avg_queries,p50_queries,p95_queries,p99_queries,early_stop_reason,"
    "query_limit_hit,query_limit_hit_frames,resource_limit_hit_frames,proven_frames,truncated_frames,"
    "stop_bound_frames,stop_exhausted_frames,run_stop_reason,runtime_seconds,runtime_kind\n";

static const char *summary_stop(const SimulationResult *result)
{
    int selected=-1;
    for (size_t i=0; i<DECODE_STOP_REASON_COUNT; ++i) {
        if (result->early_stop_counts[i]) {
            if (selected>=0)
                return "mixed";
            selected=(int)i;
        }
    }
    return selected<0 ? "none" : stop_names[selected];
}

static int write_row(FILE *stream, const DecodeConfig *decoder,
                      const SimulationConfig *config, const SimulationResult *result, bool json)
{
    double bits=(double)result->frames*result->n;
    size_t delta=strcmp(decoder->decoder_name,"lcgcd") ? 0 : decoder->delta;
    const uint64_t *stops=result->early_stop_counts;
    const char *hit=json ? (stops[DECODE_STOP_BUDGET] ? "true" : "false") : (stops[DECODE_STOP_BUDGET] ? "1" : "0");
    char runtime[64];
    if (result->runtime_available)
        snprintf(runtime,sizeof(runtime),"%.17g",result->runtime_seconds);
    else
        snprintf(runtime,sizeof(runtime),"%s",json ? "null" : "NA");

    int written=fprintf(stream,json ?
        "{\n  \"schema_version\":1,\"decoder\":\"%s\",\"delta\":%zu,\"ebn0_db\":%.17g,\"seed\":%" PRIu64
        ",\n  \"n\":%zu,\"k\":%zu,\"rate\":%.17g,\"sigma2\":%.17g,\"max_frames\":%" PRIu64
        ",\n  \"max_frame_errors\":%" PRIu64 ",\"candidate_limit\":%" PRIu64 ",\"frames\":%" PRIu64
        ",\n  \"frame_errors\":%" PRIu64 ",\"bit_errors\":%" PRIu64 ",\"hd_bit_errors\":%" PRIu64
        ",\n  \"fer\":%.17g,\"ber\":%.17g,\"hd_ber\":%.17g,\"partial_patterns_generated\":%" PRIu64
        ",\n  \"full_candidates_evaluated\":%" PRIu64 ",\"query_batches\":%" PRIu64 ",\"trellis_states_visited\":%" PRIu64
        ",\n  \"patterns_popped\":%" PRIu64 ",\"nodes_expanded\":%" PRIu64 ",\"avg_queries\":%.17g"
        ",\n  \"p50_queries\":%" PRIu64 ",\"p95_queries\":%" PRIu64 ",\"p99_queries\":%" PRIu64
        ",\n  \"early_stop_reason\":\"%s\",\"query_limit_hit\":%s,\"query_limit_hit_frames\":%" PRIu64
        ",\n  \"resource_limit_hit_frames\":%" PRIu64 ",\"proven_frames\":%" PRIu64 ",\"truncated_frames\":%" PRIu64
        ",\n  \"stop_bound_frames\":%" PRIu64 ",\"stop_exhausted_frames\":%" PRIu64
        ",\n  \"run_stop_reason\":\"%s\",\"runtime_seconds\":%s,\"runtime_kind\":\"process_cpu_decoder_calls\"\n}" :
        "1,%s,%zu,%.17g,%" PRIu64 ",%zu,%zu,%.17g,%.17g,%" PRIu64 ",%" PRIu64 ",%" PRIu64
        ",%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%.17g,%.17g,%.17g"
        ",%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%.17g"
        ",%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%s,%s,%" PRIu64 ",%" PRIu64
        ",%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%s,%s,process_cpu_decoder_calls",
        decoder->decoder_name,delta,config->ebn0_db,config->seed,result->n,result->k,(double)result->k/result->n,
        result->sigma2,config->max_frames,config->max_frame_errors,decoder->candidate_limit,
        result->frames, result->frame_errors, result->bit_errors, result->hd_bit_errors,
        (double)result->frame_errors/result->frames,result->bit_errors/bits,result->hd_bit_errors/bits,
        result->totals.partial_patterns_generated,result->totals.candidates_evaluated,result->totals.batches,
        result->totals.trellis_states_visited,result->totals.candidates_popped,result->totals.nodes_expanded,
        (double)result->totals.candidates_evaluated/result->frames,result->p50_queries,result->p95_queries,result->p99_queries,
        summary_stop(result),hit,stops[DECODE_STOP_BUDGET],stops[DECODE_STOP_RESOURCE],result->proven_frames,
        result->truncated_frames,stops[DECODE_STOP_BOUND],stops[DECODE_STOP_EXHAUSTED],
        result->error_limit_hit ? "error_limit" : "frame_limit",runtime);
    return written>=0 && !ferror(stream);
}

static int write_results(FILE *stream, const DecodeConfig *decoder, const SimulationConfig *config,
                          const SimulationResult *results, size_t count, bool json)
{
    if (!stream || !decoder || !decoder_is_known(decoder->decoder_name) || !config || !results || !results[0].frames)
        return 0;
    if (count==2 && results[0].frames!=results[1].frames)
        return 0;
    fputs(json ? "[\n" : csv_header,stream);
    for (size_t d=0; d<count; ++d) {
        DecodeConfig selected=*decoder;
        if (count==2)
            selected.decoder_name=d ? "lcgcd" : "gcd";
        if (!write_row(stream,&selected,config,&results[d],json))
            return 0;
        fputs(json && d+1<count ? ",\n" : "\n",stream);
    }
    if (json)
        fputs("]\n",stream);
    return !ferror(stream);
}

int simulation_write_csv(FILE *stream, const DecodeConfig *decoder,
                         const SimulationConfig *config, const SimulationResult *result)
{
    return write_results(stream,decoder,config,result,1,false);
}

int simulation_compare_write_csv(FILE *stream, const DecodeConfig *lc_config,
                                 const SimulationConfig *config, const SimulationResult results[2])
{
    return write_results(stream,lc_config,config,results,2,false);
}

int simulation_write_json(FILE *stream, const DecodeConfig *decoder,
                          const SimulationConfig *config, const SimulationResult *result)
{
    return write_results(stream,decoder,config,result,1,true);
}

int simulation_compare_write_json(FILE *stream, const DecodeConfig *lc_config,
                                  const SimulationConfig *config, const SimulationResult results[2])
{
    return write_results(stream,lc_config,config,results,2,true);
}
