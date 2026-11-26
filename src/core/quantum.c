#include "bplcgcd/quantum.h"

#include <stdlib.h>
#include <string.h>

/* Circular header and tail insertion adapted from cpicic_decoder
 * PrtChkMtrxAllocate/EntryInsert, parity_check_matrix.c.
 * Copyright (C) ICAIS Lab at NJU. 2020 - 2021. All rights reserved
 * Original author: Yangcan Zhou, yczhou@smail.nju.edu.cn
 * Adapted for sorted dense input, Pauli edges and offline construction.
 * See docs/third_party_sources.md for the pinned source. */
static QuantumStatus build_graph(QuantumCode *code, size_t edge_count)
{
    QbpGraph *g = &code->graph;
    g->rows = code->hx.rows;
    g->cols = code->hx.cols;
    g->edge_count = edge_count;
    if (g->rows > SIZE_MAX / sizeof(QbpEdge) ||g->cols > SIZE_MAX / sizeof(QbpEdge) ||edge_count > SIZE_MAX / sizeof(QbpEdge)) 
        return QUANTUM_INVALID_ARGUMENT;
    
    g->row_header = calloc(g->rows, sizeof(QbpEdge));
    g->col_header = calloc(g->cols, sizeof(QbpEdge));
    g->edges = edge_count ? calloc(edge_count, sizeof(QbpEdge)) : NULL;

    if (!g->row_header || !g->col_header || (edge_count && !g->edges))
        return QUANTUM_ALLOCATION_FAILED;

    //初始化行头和列头，将它们都指向自己
    for (size_t r = 0; r < g->rows; ++r) {
        QbpEdge *h = &g->row_header[r];
        *h = (QbpEdge){.row=r, .col=SIZE_MAX, .edge_index=SIZE_MAX,
                       .left=h, .right=h, .up=h, .down=h};
    }
    for (size_t c = 0; c < g->cols; ++c) {
        QbpEdge *h = &g->col_header[c];
        *h = (QbpEdge){.row=SIZE_MAX, .col=c, .edge_index=SIZE_MAX,
                       .left=h, .right=h, .up=h, .down=h};
    }
    const Pauli labels[] = {PAULI_I, PAULI_Z, PAULI_X, PAULI_Y};
    size_t k = 0;
    for (size_t r = 0; r < g->rows; ++r) {
        for (size_t c = 0; c < g->cols; ++c) {
            size_t pos = r * g->cols + c;
            Pauli p = labels[2 * code->hx.data[pos] + code->hz.data[pos]];
            if (p == PAULI_I) continue;
            QbpEdge *e = &g->edges[k];
            QbpEdge *rh = &g->row_header[r], *ch = &g->col_header[c];
            *e = (QbpEdge){.row=r, .col=c, .edge_index=k++, .pauli=p,
                           .left=rh->left, .right=rh, .up=ch->up, .down=ch}; //构建双向链表
            e->left->right = e;
            rh->left = e;
            e->up->down = e;
            ch->up = e;
        }
    }
    return QUANTUM_OK;
}

int pauli_anticommutes(Pauli a, Pauli b)
{
    return a != PAULI_I && b != PAULI_I && a != b;
}

void quantum_code_free(QuantumCode *code)
{
    if (!code) return;
    binary_matrix_free(&code->hx);
    binary_matrix_free(&code->hz);
    binary_matrix_free(&code->hsyn);
    free(code->graph.row_header);
    free(code->graph.col_header);
    free(code->graph.edges);
    *code = (QuantumCode){0};
}

QuantumStatus quantum_code_init(QuantumCode *code,
                                const BinaryMatrix *hx, const BinaryMatrix *hz)
{
    if (!code || code->hx.data || code->hz.data || code->hsyn.data ||code->graph.row_header || code->graph.col_header || code->graph.edges ||!hx || !hz || !hx->data || !hz->data || !hx->rows || !hx->cols)
        return QUANTUM_INVALID_ARGUMENT;
    if (hx->rows != hz->rows || hx->cols != hz->cols)
        return QUANTUM_DIMENSION_MISMATCH;
    
    size_t m = hx->rows, n = hx->cols;

    if (n > SIZE_MAX / 2 || m > SIZE_MAX / (2 * n)) 
        return QUANTUM_INVALID_ARGUMENT;
    size_t edge_count = 0;
    for (size_t k = 0; k < m * n; ++k) {
        if (hx->data[k] > 1 || hz->data[k] > 1) 
            return QUANTUM_INVALID_ARGUMENT;
        edge_count += (hx->data[k] | hz->data[k]) != 0;
    }
    /* Required symplectic commutation validation, offline only. */
    for (size_t a = 0; a < m; ++a) {
        for (size_t b = a + 1; b < m; ++b) {
            unsigned parity = 0;
            for (size_t c = 0; c < n; ++c)
                parity ^= (hx->data[a*n+c] & hz->data[b*n+c]) ^ (hz->data[a*n+c] & hx->data[b*n+c]);
            if (parity) 
                return QUANTUM_NONCOMMUTING;
        }
    }
    QuantumCode built = {0};
    if (binary_matrix_init(&built.hx, m, n) != GF2_OK ||binary_matrix_init(&built.hz, m, n) != GF2_OK ||binary_matrix_init(&built.hsyn, m, 2*n) != GF2_OK) {
        quantum_code_free(&built);
        return QUANTUM_ALLOCATION_FAILED;
    }
    memcpy(built.hx.data, hx->data, m*n);
    memcpy(built.hz.data, hz->data, m*n);
    for (size_t r = 0; r < m; ++r) {
        memcpy(built.hsyn.data + r*2*n, hz->data + r*n, n);
        memcpy(built.hsyn.data + r*2*n + n, hx->data + r*n, n);
    }
    QuantumStatus status = build_graph(&built, edge_count);
    if (status != QUANTUM_OK) {
        quantum_code_free(&built);
        return status;
    }
    *code = built;
    return QUANTUM_OK;
}

QuantumStatus quantum_code_syndrome(const QuantumCode *code, const Pauli *error,
                                    size_t length, BitVector *syndrome)
{
    if (!code || !code->hx.data || !error || !syndrome || !syndrome->data)
        return QUANTUM_INVALID_ARGUMENT;
    
    const QbpGraph *g = &code->graph;

    if (length != g->cols || syndrome->length != g->rows)
        return QUANTUM_DIMENSION_MISMATCH;
    for (size_t c = 0; c < length; ++c)
        if ((unsigned)error[c] > PAULI_Z) 
            return QUANTUM_INVALID_ARGUMENT;
    
    for (size_t r = 0; r < g->rows; ++r) {
        const QbpEdge *h = &g->row_header[r];
        unsigned parity = 0;
        for (const QbpEdge *e = h->right; e != h; e = e->right)
            parity ^= pauli_anticommutes(error[e->col], e->pauli);
        syndrome->data[r] = (uint8_t)parity;
    }
    return QUANTUM_OK;
}
