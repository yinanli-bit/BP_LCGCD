#include "bplcgcd/quantum_decoder.h"

#include <limits.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static DecodeStatus from_quantum(QuantumStatus s)
{
    return s == QUANTUM_OK ? DECODE_SUCCESS :
           s == QUANTUM_ALLOCATION_FAILED ? DECODE_ALLOCATION_FAILED :
           s == QUANTUM_NUMERIC_ERROR ? DECODE_NUMERIC_ERROR :
           s == QUANTUM_DIMENSION_MISMATCH ? DECODE_DIMENSION_MISMATCH : DECODE_INVALID_ARGUMENT;
}

static DecodeStatus from_gf2(Gf2Status s)
{
    return s == GF2_OK ? DECODE_SUCCESS :
           s == GF2_INCONSISTENT_SYNDROME ? DECODE_INCONSISTENT_SYNDROME :
           s == GF2_ALLOCATION_FAILED ? DECODE_ALLOCATION_FAILED : DECODE_INVALID_ARGUMENT;
}

static DecodeStatus from_slva(SlvaStatus s)
{
    return s == SLVA_OK ? DECODE_SUCCESS :
           s == SLVA_RESOURCE_LIMIT || s == SLVA_PATH_LIMIT ? DECODE_TRUNCATED :
           s == SLVA_ALLOCATION_FAILED ? DECODE_ALLOCATION_FAILED :
           s == SLVA_NUMERIC_ERROR ? DECODE_NUMERIC_ERROR : DECODE_INVALID_ARGUMENT;
}

void quantum_decoder_free(QuantumDecoderWorkspace *w)
{
    if (!w) return;
    slva_free(&w->list);
    qbp_workspace_free(&w->bp);
    binary_systematic_free(&w->syndrome_form);
    binary_systematic_free(&w->stabilizer_form);
    BitVector *vectors[] = {&w->reference,&w->residual,&w->transformed,&w->system,
                            &w->partial,&w->binary,&w->syndrome,&w->local_rhs};
    for (size_t i = 0; i < sizeof(vectors)/sizeof(*vectors); ++i) bit_vector_free(vectors[i]);
    void *arrays[] = {w->swt,w->partial_swt,w->reliability,w->candidate,w->recovery,
                     w->local_p,w->local.basis_order,w->trellis.column_masks,w->trellis.original_columns};
    for (size_t i = 0; i < sizeof(arrays)/sizeof(*arrays); ++i) free(arrays[i]);
    *w = (QuantumDecoderWorkspace){0};
}

DecodeStatus quantum_decoder_init(QuantumDecoderWorkspace *w,
                                  const QuantumCode *code, const QuantumLcgcdConfig *config)
{
    if (!w || w->code || !code || !code->hx.data || !config ||(config->post_policy != BP_POST_ON_FAILURE && config->post_policy != BP_POST_ALWAYS) ||(config->stop_policy != QUANTUM_STOP_TRIVIAL && config->stop_policy != QUANTUM_STOP_DISABLED && config->stop_policy != QUANTUM_STOP_DAI))
        return DECODE_INVALID_ARGUMENT;
    
    size_t n = code->graph.cols, m = code->graph.rows, N = code->hsyn.cols;

    if (n > SIZE_MAX/sizeof(PauliSwt) || N > SIZE_MAX/sizeof(double) || N > SIZE_MAX/sizeof(size_t)) 
        return DECODE_ALLOCATION_FAILED;
    
    QuantumDecoderWorkspace built = {.code=code, .config=*config};
    DecodeStatus status = from_gf2(binary_matrix_systematize(&code->hsyn, &built.syndrome_form));
    if (status != DECODE_SUCCESS) goto fail;
    /* S and Hsyn have different column meanings. Only S defines equivalence. */
    BinaryMatrix stabilizer = {0};
    status = from_gf2(binary_matrix_init(&stabilizer, m, N));
    if (status != DECODE_SUCCESS) goto fail;
    for (size_t r = 0; r < m; ++r) {
        memcpy(stabilizer.data+r*N, code->hx.data+r*n, n);
        memcpy(stabilizer.data+r*N+n, code->hz.data+r*n, n);
    }
    status = from_gf2(binary_matrix_systematize(&stabilizer, &built.stabilizer_form));
    binary_matrix_free(&stabilizer);
    if (status != DECODE_SUCCESS) goto fail;
    size_t rank = built.syndrome_form.rank, delta = config->delta, k = N-rank;
    if (delta > rank) {
         status=DECODE_DIMENSION_MISMATCH; goto fail; 
    }
    if (delta >= sizeof(size_t)*CHAR_BIT || ((size_t)1 << delta) > (config->state_limit ? config->state_limit : 4096)) {
        status = DECODE_TRUNCATED;
        goto fail;
    }
    size_t stages = k+delta;
    built.swt = malloc(n*sizeof(PauliSwt));
    built.partial_swt = malloc(n*sizeof(PartialSwt));
    built.reliability = malloc(N*sizeof(double));
    built.candidate = malloc(n*sizeof(Pauli));
    built.recovery = malloc(n*sizeof(Pauli));
    built.local_p = rank && k ? malloc(rank*k) : NULL;
    built.local.basis_order = rank ? malloc(rank*sizeof(size_t)) : NULL;
    built.trellis = (SyndromeTrellis){.length=stages, .column_masks=malloc(stages*sizeof(size_t)), .original_columns=malloc(stages*sizeof(size_t))};

    status = DECODE_ALLOCATION_FAILED;
    if (!built.swt || !built.partial_swt || !built.reliability || !built.candidate || !built.recovery ||(rank && k && !built.local_p) || (rank && !built.local.basis_order) || !built.trellis.column_masks || !built.trellis.original_columns) 
        goto fail;
    BitVector *vectors[] = {&built.reference,&built.residual,&built.transformed,&built.system,
                            &built.partial,&built.binary,&built.syndrome,&built.local_rhs};
    const size_t lengths[] = {N,m,m,N,stages,N,m,rank};
    for (size_t i = 0; i < sizeof(vectors)/sizeof(*vectors); ++i)
        if (lengths[i] && bit_vector_init(vectors[i], lengths[i]) != GF2_OK) 
            goto fail;
    built.local.p_a = (BinaryMatrix){delta,k,built.local_p};
    built.local.p_ac = (BinaryMatrix){rank-delta,k,built.local_p ? built.local_p+delta*k : NULL};
    built.local.s_a = (BitVector){delta,built.local_rhs.data};
    built.local.s_ac = (BitVector){rank-delta,built.local_rhs.data ? built.local_rhs.data+delta : NULL};
    SlvaLimits limits = {.max_candidates=config->slva_max_candidates ? config->slva_max_candidates : 65536,
        .max_bytes=config->slva_max_bytes ? config->slva_max_bytes : 64*1024*1024};
    status = from_slva(slva_workspace_init(&built.list, stages, (size_t)1 << delta, &limits));
    if (status != DECODE_SUCCESS) goto fail;
    status = from_quantum(qbp_workspace_init(&built.bp, code));
    if (status != DECODE_SUCCESS) goto fail;
    *w = built;
    return DECODE_SUCCESS;
fail:
    quantum_decoder_free(&built);
    return status;
}


static DecodeStatus begin_frame(QuantumDecoderWorkspace *w, const BitVector *measured, QuantumDecodeResult *out)
{
    if (!out) 
        return DECODE_INVALID_ARGUMENT;
    
    *out = (QuantumDecodeResult){0};

    if (!w || !w->code || !measured || !measured->data) 
        return out->status=DECODE_INVALID_ARGUMENT;
    if (measured->length != w->code->graph.rows) 
        return out->status=DECODE_DIMENSION_MISMATCH;
    for (size_t r = 0; r < measured->length; ++r)
        if (measured->data[r] > 1) 
            return out->status=DECODE_INVALID_ARGUMENT;
    
    out->length = w->code->graph.cols;
    out->recovery = w->recovery;
    return DECODE_SUCCESS;
}

static DecodeStatus prepare_posterior(QuantumDecoderWorkspace *w, const QbpResult *post,
                                      QuantumDecodeResult *out)
{
    if (!post) return DECODE_INVALID_ARGUMENT;
    if (post->length != out->length) return DECODE_DIMENSION_MISMATCH;
    out->posterior = *post;
    return from_quantum(qbp_posterior_to_swt(post, w->swt));
}

static double pauli_cost(const QuantumDecoderWorkspace *w, const Pauli *candidate)
{
    double score = 0;
    for (size_t i = 0; i < w->code->graph.cols; ++i) score += w->swt[i].w[candidate[i]];
    return score;
}

/*add_up()、sum_error()、nextafterl()：为了确保浮点误差不会导致错误提前停止。*/
static long double add_up(long double a, long double b)
{
    return nextafterl(a+b, INFINITY); 
}

/* Nonnegative double sums: |fl(sum)-sum| <= gamma_count*sum.
 * Upward evaluation also covers platforms where long double == double. */
static long double sum_error(size_t count, long double total)
{
    long double nu = nextafterl((long double)count*(DBL_EPSILON/2), INFINITY);
    if (nu >= 1) return INFINITY;
    long double gamma = nextafterl(nu/nextafterl(1-nu, -INFINITY), INFINITY);
    return nextafterl(gamma*total, INFINITY); //把求和结果朝 \(+\infty\) 推一点
}

/* 2025 Eq. (23) needs a Pauli lower bound, not raw marginal cost.
 * For each qubit C_i=max_P(B_i(P)-w_i(P)), so W(Q)>=B_M(Q)-sum C_i.
 * See docs/quantum_early_stopping.md for the proof and fixed-basis DAI mapping. */
static void prepare_stopping(QuantumDecoderWorkspace *w)
{
    size_t n = w->code->graph.cols, N = 2*n;
    memset(w->binary.data, 0, N);
    long double binary_total = 0, pauli_total = 0, correction = 0;
    for (size_t j = 0; j < w->trellis.length; ++j) {
        size_t col = w->trellis.original_columns[j];
        w->binary.data[col] = 1;
        binary_total = add_up(binary_total, w->reliability[col]);
    }
    static const uint8_t x[] = {0,1,1,0}, z[] = {0,0,1,1};
    w->dai_threshold = 0;
    for (size_t i = 0; i < n; ++i) {
        long double worst = 0;
        double largest = 0, min1 = INFINITY, min2 = INFINITY;
        for (size_t p = 0; p < 4; ++p) {
            long double bx = w->binary.data[i] && (x[p] != w->reference.data[i]) ? w->reliability[i] : 0;
            long double bz = w->binary.data[n+i] && (z[p] != w->reference.data[n+i]) ? w->reliability[n+i] : 0;
            long double gap = nextafterl(add_up(bx,bz)-w->swt[i].w[p], INFINITY);
            worst = fmaxl(worst, gap);
            double weight = w->swt[i].w[p];
            largest = fmax(largest, weight);
            if (weight < min1) { min2=min1; min1=weight; }
            else if (weight < min2) min2=weight;
        }
        correction = add_up(correction, worst);
        pauli_total = add_up(pauli_total, largest);
        /* LRB is the unique qubit projection of Ac; count a mixed/Y qubit once. */
        if (w->config.stop_policy == QUANTUM_STOP_DAI &&
            (!w->binary.data[i] || !w->binary.data[n+i])) {
            double r = sqrt(min2), e = exp(-r); /* Eq. (26): second-smallest sqrt(weight). */
            w->dai_threshold += r*e/(1+e);       /* Eq. (25), overflow-safe. */
        }
    }
    w->stop_correction = add_up(correction, add_up(sum_error(w->trellis.length,binary_total),
                                                  sum_error(n,pauli_total)));
}

static double remaining_bound(const QuantumDecoderWorkspace *w, double binary_next)
{
    long double lower = nextafterl((long double)binary_next-w->stop_correction, -INFINITY);
    if (lower <= 0) return 0;
    double bound = (double)lower;
    return (long double)bound > lower ? nextafter(bound,-INFINITY) : bound;
}

static bool stop_search(QuantumDecoderWorkspace *w, double bound, QuantumDecodeResult *out)
{
    out->remaining_lower_bound = bound;
    if (bound >= out->cost) {
        out->optimality_proven = true;
        out->early_stop = QUANTUM_EARLY_TRIVIAL;
        out->stop_reason = DECODE_STOP_BOUND;
        return true;
    }
    if (w->config.stop_policy == QUANTUM_STOP_DAI && bound+out->dai_threshold >= out->cost) {
        out->early_stop = QUANTUM_EARLY_DAI;
        /* The positive DAI threshold estimates omitted cost; it is not a bound. */
        return true;
    }
    return false;
}

static DecodeStatus postprocess(QuantumDecoderWorkspace *w, const BitVector *measured,
                                QuantumDecodeResult *out)
{
    out->postprocess_invoked = true;
    DecodeStatus status = from_quantum(quantum_swt_to_partial(w->swt, out->length, w->partial_swt, &w->reference, w->reliability));

    if (status != DECODE_SUCCESS) 
        return out->status=status;
    status = from_gf2(binary_matrix_multiply_vector(&w->code->hsyn, &w->reference, &w->residual));
    if (status != DECODE_SUCCESS) 
        return out->status=status;
    for (size_t r = 0; r < measured->length; ++r) 
        w->residual.data[r] ^= measured->data[r];
    
    const BinarySystematic *form = &w->syndrome_form;
    size_t N = form->reduced.cols, rank = form->rank, k = N-rank, delta = w->config.delta;
    status = from_gf2(local_constraints_prepare(form, w->reliability, N, &w->residual,
                                               delta, &w->local, &w->transformed));
    if (status != DECODE_SUCCESS) 
        return out->status=status;
    TrellisStatus prepared = syndrome_trellis_prepare(form, &w->local, w->reliability, N,w->config.state_limit ? w->config.state_limit : 4096, &w->trellis);

    if (prepared != TRELLIS_OK)
        return out->status=prepared == TRELLIS_NUMERIC_ERROR ? DECODE_NUMERIC_ERROR : DECODE_INVALID_ARGUMENT;
    bool stopping = w->config.stop_policy != QUANTUM_STOP_DISABLED;
    if (stopping) {
        prepare_stopping(w);
        out->dai_threshold = w->dai_threshold;
    }
    status = from_slva(slva_workspace_reset(&w->list, &w->trellis));
    if (status != DECODE_SUCCESS) return out->status=status;
    static const Pauli labels[] = {PAULI_I,PAULI_Z,PAULI_X,PAULI_Y};
    for (;;) {
        /* Nonnegative Pauli SWT proves a zero-score incumbent immediately,
         * without expanding SLVA or consuming another candidate/resource slot. */
        if (stopping && out->candidate_valid && out->cost == 0 && stop_search(w,0,out)) {
            status = DECODE_SUCCESS;
            break;
        }
        double ordering_cost;
        SlvaStatus next = slva_peek(&w->list, &ordering_cost);
        if (next == SLVA_DONE) {
            out->optimality_proven = true;
            out->stop_reason = DECODE_STOP_EXHAUSTED;
            status = DECODE_SUCCESS;
            break;
        }
        if (next != SLVA_OK) {
            status = from_slva(next);
            if (status == DECODE_TRUNCATED) out->stop_reason = DECODE_STOP_RESOURCE;
            break;
        }
        if (stopping && out->candidate_valid && stop_search(w,remaining_bound(w,ordering_cost),out)) {
            status = DECODE_SUCCESS;
            break;
        }
        if (w->config.candidate_limit && out->stats.candidates_evaluated >= w->config.candidate_limit) {
            status = DECODE_TRUNCATED;
            out->stop_reason = DECODE_STOP_BUDGET;
            break;
        }
        next = slva_next(&w->list, &w->partial, &ordering_cost);
        if (next != SLVA_OK) { status=from_slva(next); break; }
        ++out->stats.candidates_popped;
        /* LC-GCD Eq. (9)/(10), same fixed-basis reconstruction as classical. */
        memset(w->system.data, 0, rank);
        memcpy(w->system.data+rank, w->partial.data+delta, k);
        for (size_t i = 0; i < delta; ++i)
            w->system.data[form->orig_to_sys[w->local.basis_order[i]]] = w->partial.data[i];
        for (size_t i = 0; i < rank-delta; ++i) {
            unsigned bit = w->local.s_ac.data[i];
            for (size_t j = 0; j < k; ++j) 
                bit ^= w->local.p_ac.data[i*k+j] & w->system.data[rank+j];
            w->system.data[form->orig_to_sys[w->local.basis_order[delta+i]]] = (uint8_t)bit;
        }
        for (size_t j = 0; j < N; ++j)
            w->binary.data[form->sys_to_orig[j]] = w->system.data[j] ^ w->reference.data[form->sys_to_orig[j]]; //N位，[E_x,E_z]
        for (size_t i = 0; i < out->length; ++i)
            w->candidate[i] = labels[2*w->binary.data[i]+w->binary.data[out->length+i]]; //对应转化成n位puali错误模式
        /* Raw ordering cost never replaces complete Pauli scoring. */
        double cost = pauli_cost(w, w->candidate); //完整的Pauli评分
        if (!isfinite(cost)) { status=DECODE_NUMERIC_ERROR; break; }
        ++out->stats.candidates_evaluated;
        if (!out->candidate_valid || cost < out->cost) {
            memcpy(w->recovery, w->candidate, out->length*sizeof(Pauli));
            out->cost = cost;
            out->candidate_valid = true;
        }
    }
    out->stats.nodes_expanded = w->list.nodes_expanded;
    out->stats.partial_patterns_generated = w->list.next_sequence;

    for (size_t layer = 0; layer < w->trellis.length; ++layer)
        for (size_t state = 0; state < w->trellis.state_count; ++state)
            out->stats.trellis_states_visited += isfinite(w->list.cells[layer*w->trellis.state_count+state].metric);
    if (out->candidate_valid) {
        quantum_code_syndrome(w->code, w->recovery, out->length, &w->syndrome);
        out->syndrome_matched = memcmp(w->syndrome.data, measured->data, measured->length) == 0;
        if (!out->syndrome_matched) status = DECODE_INVALID_CANDIDATE;
    }
    if (status != DECODE_SUCCESS && status != DECODE_TRUNCATED)
        out->candidate_valid = out->syndrome_matched = out->optimality_proven = false;
    return out->status=status;
}

DecodeStatus quantum_lcgcd_decode(QuantumDecoderWorkspace *w,
                                  const QbpResult *posterior, const BitVector *measured,
                                  QuantumDecodeResult *out)
{
    DecodeStatus status = begin_frame(w, measured, out);
    if (status != DECODE_SUCCESS) return status;
    status = prepare_posterior(w, posterior, out);
    return status == DECODE_SUCCESS ? postprocess(w, measured, out) : (out->status=status);
}

DecodeStatus quantum_bp_lcgcd_decode(QuantumDecoderWorkspace *w,
                                     const QbpConfig *bp_config, const QbpPrior *prior,
                                     const BitVector *measured, QuantumDecodeResult *out)
{
    DecodeStatus status = begin_frame(w, measured, out);
    if (status != DECODE_SUCCESS) return status;

    QbpResult posterior;
    status = from_quantum(qbp_decode(&w->bp, bp_config, prior, measured, &posterior));

    if (status != DECODE_SUCCESS) return out->status=status;

    status = prepare_posterior(w, &posterior, out);

    if (status != DECODE_SUCCESS) return out->status=status;
    if (posterior.syndrome_matched && w->config.post_policy == BP_POST_ON_FAILURE) {
        memcpy(w->recovery, posterior.hard_estimate, out->length*sizeof(Pauli));
        out->cost = pauli_cost(w, w->recovery);
        if (!isfinite(out->cost)) return out->status=DECODE_NUMERIC_ERROR;
        out->candidate_valid = out->syndrome_matched = true;
        return out->status=DECODE_SUCCESS;
    }
    return postprocess(w, measured, out);
}

QuantumStatus quantum_recovery_check(QuantumDecoderWorkspace *w,
                                     const Pauli *physical, const Pauli *recovery,
                                     size_t length, QuantumRecoveryCheck *out)
{
    if (!out) return QUANTUM_INVALID_ARGUMENT;

    *out = (QuantumRecoveryCheck){0};

    if (!w || !w->code || !physical || !recovery) 
        return QUANTUM_INVALID_ARGUMENT;
    if (length != w->code->graph.cols) 
        return QUANTUM_DIMENSION_MISMATCH;
    
    static const uint8_t x[] = {0,1,1,0}, z[] = {0,0,1,1};
    static const Pauli labels[] = {PAULI_I,PAULI_Z,PAULI_X,PAULI_Y};
    for (size_t i = 0; i < length; ++i) {
        if ((unsigned)physical[i] > PAULI_Z || (unsigned)recovery[i] > PAULI_Z)
            return QUANTUM_INVALID_ARGUMENT;
        w->binary.data[i] = x[physical[i]] ^ x[recovery[i]];
        w->binary.data[length+i] = z[physical[i]] ^ z[recovery[i]];
        w->candidate[i] = labels[2*w->binary.data[i]+w->binary.data[length+i]];
    }
    //第一步：先确保recovery修对了measured syndrome
    quantum_code_syndrome(w->code, w->candidate, length, &w->syndrome);
    out->syndrome_matched = true;
    for (size_t r = 0; r < w->syndrome.length; ++r)
        if (w->syndrome.data[r]) {
            out->syndrome_matched=false; 
            return QUANTUM_OK; 
        }
    /* D in rowspan(S): after the OFFLINE column permutation, D_free=D_basis*P(系统形矩阵的行空间性质) */
    //syndrome一致不一定译对，第二步判断真实error和译码error pattern的残差是不是一个stabilizer
    const BinarySystematic *form = &w->stabilizer_form;//[H_x|H_z]
    size_t N = form->reduced.cols, rank = form->rank;
    for (size_t j = rank; j < N; ++j) {
        unsigned bit = w->binary.data[form->sys_to_orig[j]];
        for (size_t i = 0; i < rank; ++i)
            bit ^= w->binary.data[form->sys_to_orig[i]] & form->reduced.data[i*N+j];
        if (bit) return QUANTUM_OK;
    }
    out->logical_success = true;
    return QUANTUM_OK;
}
