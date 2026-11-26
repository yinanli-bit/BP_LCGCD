#include "bplcgcd/qbp.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

void qbp_workspace_free(QbpWorkspace *w)
{
    if (!w) return;
    free(w->messages);
    free(w->hard_estimate);
    free(w->syndrome);
    *w = (QbpWorkspace){0};
}

QuantumStatus qbp_workspace_init(QbpWorkspace *w, const QuantumCode *code)
{
    if (!w || w->code || w->messages || w->hard_estimate || w->syndrome ||
        !code || !code->hx.data) return QUANTUM_INVALID_ARGUMENT;
    size_t n = code->graph.cols, m = code->graph.rows, edges = code->graph.edge_count;
    size_t capacity = SIZE_MAX / sizeof(double);
    if (n > capacity / 3 || edges > (capacity - 3*n) / 5 ||
        n > SIZE_MAX / sizeof(Pauli)) return QUANTUM_INVALID_ARGUMENT;
    QbpWorkspace built = {.code=code};
    built.messages = calloc(5*edges + 3*n, sizeof(double));
    built.hard_estimate = malloc(n * sizeof(Pauli));
    built.syndrome = malloc(m);
    if (!built.messages || !built.hard_estimate || !built.syndrome) {
        qbp_workspace_free(&built);
        return QUANTUM_ALLOCATION_FAILED;
    }
    built.v2c_x = built.messages;
    built.v2c_y = built.v2c_x + edges;
    built.v2c_z = built.v2c_y + edges;
    built.c2v = built.v2c_z + edges;
    built.commutation = built.c2v + edges;
    built.post_x = built.commutation + edges;
    built.post_y = built.post_x + n;
    built.post_z = built.post_y + n;
    *w = built;
    return QUANTUM_OK;
}

static double log_add_exp(double a, double b)
{
    double high = fmax(a, b), low = fmin(a, b);
    return high + log1p(exp(low - high));//采用数值稳定写法，避免大正数 exp() 溢出，实际计算的是log(exp(a) + exp(b))
}

/* 2025 Eq. (6): commuting mass {I,P} versus anticommuting mass {Q,R}.
 * Choosing indexes here keeps the reduction out of check-node arithmetic. */
static double commutation_llr(const double xyz[3], Pauli p)
{
    size_t own = (size_t)p - 1;
    size_t q = (own + 1) % 3, r = (own + 2) % 3;
    return log_add_exp(0, -xyz[own]) - log_add_exp(-xyz[q], -xyz[r]);
}

double qbp_commutation_llr(const double msg_xyz[3], Pauli stabilizer_pauli)
{
    if (!msg_xyz || stabilizer_pauli < PAULI_X || stabilizer_pauli > PAULI_Z ||
        !isfinite(msg_xyz[0]) || !isfinite(msg_xyz[1]) || !isfinite(msg_xyz[2]))
        return NAN;
    return commutation_llr(msg_xyz, stabilizer_pauli);
}

/* min1/min2 and two-pass sparse check processing adapted from
 * cpicic_decoder/dec.c: CheckNodeProcMinSum, Min1Min2Scaled.
 * Copyright (C) ICAIS Lab at NJU. 2020 - 2021. All rights reserved
 * Original author: Yangcan Zhou, yczhou@smail.nju.edu.cn
 * Adapted: scalar commutation LLR from XYZ, syndrome sign, edge-indexed
 * workspace, floating clip. Quantum formulas: 2025 Eq. (6), (12).
 * Source commit and differences: docs/third_party_sources.md. */
static QuantumStatus check_update(QbpWorkspace *w, const BitVector *measured,
                                   double nms_scale, double clip)
{
    const QbpGraph *g = &w->code->graph;
    for (size_t row = 0; row < g->rows; ++row) {
        const QbpEdge *h = &g->row_header[row];
        double min1 = clip, min2 = clip;
        size_t min_index = SIZE_MAX;
        unsigned parity = measured->data[row];
        for (const QbpEdge *e = h->right; e != h; e = e->right) {
            size_t k = e->edge_index;
            const double xyz[] = {w->v2c_x[k], w->v2c_y[k], w->v2c_z[k]};
            double lambda = commutation_llr(xyz, e->pauli);
            if (!isfinite(lambda)) 
                return QUANTUM_NUMERIC_ERROR;
            w->commutation[k] = lambda;
            double magnitude = fabs(lambda);
            parity ^= lambda < 0;
            if (magnitude < min1) {
                min2 = min1;
                min1 = magnitude;
                min_index = k;
            } else if (magnitude < min2) {
                min2 = magnitude;
            }
        }
        for (const QbpEdge *e = h->right; e != h; e = e->right) {
            size_t k = e->edge_index;
            double magnitude = nms_scale * (k == min_index ? min2 : min1);
            unsigned negative = parity ^ (w->commutation[k] < 0);
            w->c2v[k] = negative ? -magnitude : magnitude;
        }
    }
    return QUANTUM_OK;
}

static QuantumStatus posterior_update(QbpWorkspace *w, const QbpPrior *prior)
{
    const QbpGraph *g = &w->code->graph;
    for (size_t c = 0; c < g->cols; ++c) {
        double x = prior->prior_x[c], y = prior->prior_y[c], z = prior->prior_z[c];
        const QbpEdge *h = &g->col_header[c];
        for (const QbpEdge *e = h->down; e != h; e = e->down) {
            double msg = w->c2v[e->edge_index];
            if (e->pauli != PAULI_X) x += msg;
            if (e->pauli != PAULI_Y) y += msg;
            if (e->pauli != PAULI_Z) z += msg;
        }
        if (!isfinite(x) || !isfinite(y) || !isfinite(z)) 
            return QUANTUM_NUMERIC_ERROR;
        w->post_x[c] = x;
        w->post_y[c] = y;
        w->post_z[c] = z;
        /* Eq. (8); ties choose first of I, X, Y, Z. */
        double best = 0;
        Pauli hard = PAULI_I;
        if (x < best) { best = x; hard = PAULI_X; }
        if (y < best) { best = y; hard = PAULI_Y; }
        if (z < best) { hard = PAULI_Z; }
        w->hard_estimate[c] = hard;
    }
    return QUANTUM_OK;
}

static QuantumStatus variable_update(QbpWorkspace *w)
{
    const QbpGraph *g = &w->code->graph;
    for (size_t c = 0; c < g->cols; ++c) {
        const QbpEdge *h = &g->col_header[c];
        for (const QbpEdge *e = h->down; e != h; e = e->down) {
            size_t k = e->edge_index;
            double msg = w->c2v[k];
            /* Eq. (10): remove only this check's anticommuting terms. */
            w->v2c_x[k] = w->post_x[c] - (e->pauli != PAULI_X ? msg : 0);
            w->v2c_y[k] = w->post_y[c] - (e->pauli != PAULI_Y ? msg : 0);
            w->v2c_z[k] = w->post_z[c] - (e->pauli != PAULI_Z ? msg : 0);
            if (!isfinite(w->v2c_x[k]) || !isfinite(w->v2c_y[k]) || !isfinite(w->v2c_z[k]))
                return QUANTUM_NUMERIC_ERROR;
        }
    }
    return QUANTUM_OK;
}

QuantumStatus qbp_decode(QbpWorkspace *w, const QbpConfig *config,
                         const QbpPrior *prior, const BitVector *measured,
                         QbpResult *result)
{
    if (!result) 
        return QUANTUM_INVALID_ARGUMENT;
    *result = (QbpResult){0};
    if (!w || !w->code || !config || !prior || !prior->prior_x ||!prior->prior_y || !prior->prior_z || !measured || !measured->data ||!config->max_iterations || !isfinite(config->nms_scale) || config->nms_scale <= 0 || config->nms_scale > 1 || !isfinite(config->message_clip) || config->message_clip < 0)
        return QUANTUM_INVALID_ARGUMENT;
    const QbpGraph *g = &w->code->graph;
    if (prior->length != g->cols || measured->length != g->rows)
        return QUANTUM_DIMENSION_MISMATCH;
    for (size_t c = 0; c < g->cols; ++c)
        if (!isfinite(prior->prior_x[c]) || !isfinite(prior->prior_y[c]) || !isfinite(prior->prior_z[c])) 
            return QUANTUM_INVALID_ARGUMENT;
    for (size_t r = 0; r < g->rows; ++r)
        if (measured->data[r] > 1) 
            return QUANTUM_INVALID_ARGUMENT;
    for (size_t k = 0; k < g->edge_count; ++k) {
        size_t c = g->edges[k].col;
        w->v2c_x[k] = prior->prior_x[c];
        w->v2c_y[k] = prior->prior_y[c];
        w->v2c_z[k] = prior->prior_z[c];
        w->c2v[k] = 0;
    }
    double clip = config->message_clip ? config->message_clip : 50;
    BitVector syndrome = {g->rows, w->syndrome};
    /* Flooding: complete all checks before any posterior/V2C updates. */
    for (size_t iteration = 0; iteration < config->max_iterations; ++iteration) {
        QuantumStatus status = check_update(w, measured, config->nms_scale, clip);
        if (status != QUANTUM_OK) 
            return status;
        
        status = posterior_update(w, prior);
        if (status != QUANTUM_OK) 
            return status;
        
        status = quantum_code_syndrome(w->code, w->hard_estimate, g->cols, &syndrome);
        if (status != QUANTUM_OK) 
            return status;
        int matched = memcmp(w->syndrome, measured->data, g->rows) == 0;
        if (matched || iteration == config->max_iterations - 1) {
            *result = (QbpResult){.length=g->cols, .iterations_used=iteration+1,
                .syndrome_matched=matched, .hard_estimate=w->hard_estimate,
                .post_x=w->post_x, .post_y=w->post_y, .post_z=w->post_z};
            return QUANTUM_OK;
        }

        status = variable_update(w);
        if (status != QUANTUM_OK) return status;
    }
    return QUANTUM_OK;
}
