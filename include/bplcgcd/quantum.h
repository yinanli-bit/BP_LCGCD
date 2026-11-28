#ifndef BPLCGCD_QUANTUM_H
#define BPLCGCD_QUANTUM_H

#include "bplcgcd/gf2.h"

typedef enum {
    PAULI_I = 0, PAULI_X, PAULI_Y, PAULI_Z
} Pauli;

typedef enum {
    QUANTUM_OK = 0,
    QUANTUM_INVALID_ARGUMENT,
    QUANTUM_DIMENSION_MISMATCH,
    QUANTUM_NONCOMMUTING,
    QUANTUM_ALLOCATION_FAILED,
    QUANTUM_NUMERIC_ERROR
} QuantumStatus;

/* Four-direction topology adapted from cpicic_decoder Entry/PrtChkMtrx.
 * Copyright (C) ICAIS Lab at NJU. 2020 - 2021. All rights reserved
 * Original author: Yangcan Zhou, yczhou@smail.nju.edu.cn
 * Adapted: Pauli labels, stable edge_index, separate message workspace.
 * Source commit and changes: docs/third_party_sources.md. */
typedef struct QbpEdge {
    size_t row, col, edge_index;
    Pauli pauli;
    struct QbpEdge *left, *right, *up, *down;
} QbpEdge;

typedef struct {
    size_t rows, cols, edge_count;
    QbpEdge *row_header, *col_header, *edges;
} QbpGraph;

/* **hx,hz,hsyn,graph**: 
 * Owning offline object; initialize with {0}, do not copy by value or mutate
 * its members after construction. Inputs are copied. Headers are circular:
 * traverse right/down until returning to the corresponding header. */
typedef struct {
    BinaryMatrix hx, hz; /* Stabilizer S=[Hx|Hz], both m by n. */
    BinaryMatrix hsyn;   /* Syndrome matrix [Hz|Hx], m by 2n. */
    QbpGraph graph;      /* One edge per non-I Pauli, including Y. */
} QuantumCode;

/* Pauli arguments must be I/X/Y/Z; phases are ignored. */
int pauli_anticommutes(Pauli a, Pauli b);
QuantumStatus quantum_code_init(QuantumCode *code,
                                const BinaryMatrix *hx, const BinaryMatrix *hz);
void quantum_code_free(QuantumCode *code);
/* Caller-owned output; original check/qubit order, no allocation. */
QuantumStatus quantum_code_syndrome(const QuantumCode *code, const Pauli *error,
                                    size_t length, BitVector *syndrome);

#endif
