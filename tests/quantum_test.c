#include "bplcgcd/quantum_channel.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

static unsigned checks;
#define CHECK(c) do { if (!(c)) { \
    fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, #c); exit(EXIT_FAILURE); \
} ++checks; } while (0)

static int near(double a, double b) { return fabs(a-b) < 1e-12; }

static void test_pauli_and_graph(void)
{
    const unsigned anti[4][4] = {{0,0,0,0}, {0,0,1,1}, {0,1,0,1}, {0,1,1,0}};
    for (unsigned a = 0; a < 4; ++a)
        for (unsigned b = 0; b < 4; ++b)
            CHECK(pauli_anticommutes((Pauli)a, (Pauli)b) == (int)anti[a][b]);

    /* Mixed, commuting generators; product row, empty row, isolated column. */
    const Pauli stabilizer[4][5] = {
        {PAULI_X,PAULI_Y,PAULI_Z,PAULI_I,PAULI_I},
        {PAULI_Z,PAULI_X,PAULI_I,PAULI_Y,PAULI_I},
        {PAULI_Y,PAULI_Z,PAULI_Z,PAULI_Y,PAULI_I},
        {PAULI_I,PAULI_I,PAULI_I,PAULI_I,PAULI_I}
    };
    uint8_t hx_data[20], hz_data[20];
    for (size_t r = 0; r < 4; ++r) {
        for (size_t c = 0; c < 5; ++c) {
            Pauli p = stabilizer[r][c];
            hx_data[r*5+c] = p == PAULI_X || p == PAULI_Y;
            hz_data[r*5+c] = p == PAULI_Z || p == PAULI_Y;
        }
    }
    BinaryMatrix hx = {4,5,hx_data}, hz = {4,5,hz_data};
    QuantumCode code = {0};
    CHECK(quantum_code_init(&code, &hx, &hz) == QUANTUM_OK);
    CHECK(code.graph.edge_count == 10);
    CHECK(code.hsyn.rows == 4 && code.hsyn.cols == 10);
    CHECK(code.hx.data != hx.data && code.hz.data != hz.data);
    for (size_t r = 0; r < 4; ++r) {
        for (size_t c = 0; c < 5; ++c) {
            CHECK(code.hsyn.data[r*10+c] == hz_data[r*5+c]);
            CHECK(code.hsyn.data[r*10+c+5] == hx_data[r*5+c]);
        }
    }
    unsigned row_seen[10] = {0}, col_seen[10] = {0};
    for (size_t r = 0; r < 4; ++r) {
        const QbpEdge *h = &code.graph.row_header[r];
        size_t forward = 0, backward = 0, last_col = 0;
        for (const QbpEdge *e = h->right; e != h; e = e->right) {
            CHECK(forward++ < 5 && e->row == r);
            CHECK(e->pauli == stabilizer[r][e->col] && e->pauli != PAULI_I);
            CHECK(e->left->right == e && e->right->left == e);
            CHECK(e->up->down == e && e->down->up == e);
            CHECK(e->edge_index < 10 && &code.graph.edges[e->edge_index] == e);
            if (forward > 1) CHECK(last_col < e->col);
            last_col = e->col;
            ++row_seen[e->edge_index];
        }
        for (const QbpEdge *e = h->left; e != h; e = e->left) CHECK(++backward <= forward);
        CHECK(forward == backward);
    }
    for (size_t c = 0; c < 5; ++c) {
        const QbpEdge *h = &code.graph.col_header[c];
        size_t forward = 0, backward = 0, last_row = 0;
        for (const QbpEdge *e = h->down; e != h; e = e->down) {
            CHECK(forward++ < 4 && e->col == c);
            CHECK(e->up->down == e && e->down->up == e);
            if (forward > 1) CHECK(last_row < e->row);
            last_row = e->row;
            ++col_seen[e->edge_index];
        }
        for (const QbpEdge *e = h->up; e != h; e = e->up) CHECK(++backward <= forward);
        CHECK(forward == backward);
    }
    for (size_t k = 0; k < 10; ++k) CHECK(row_seen[k] == 1 && col_seen[k] == 1);
    CHECK(code.graph.row_header[3].right == &code.graph.row_header[3]);
    CHECK(code.graph.col_header[4].down == &code.graph.col_header[4]);
    /* Constructor owns copies; later changes to caller matrices have no effect. */
    memset(hx_data, 0, sizeof(hx_data));
    memset(hz_data, 0, sizeof(hz_data));
    CHECK(code.hx.data[0] == 1 && code.hz.data[0] == 0);

    uint8_t syn_data[4], gf2_data[4], bits[10];
    BitVector syn = {4,syn_data}, gf2_syn = {4,gf2_data}, binary_error = {10,bits};
    Pauli error[5];
    /* Independent dense symplectic oracle over all 4^5 Pauli errors. */
    for (unsigned word = 0; word < 1024; ++word) {
        for (size_t c = 0; c < 5; ++c) {
            error[c] = (Pauli)((word >> (2*c)) & 3u);
            bits[c] = error[c] == PAULI_X || error[c] == PAULI_Y;
            bits[c+5] = error[c] == PAULI_Z || error[c] == PAULI_Y;
        }
        CHECK(quantum_code_syndrome(&code, error, 5, &syn) == QUANTUM_OK);
        CHECK(binary_matrix_multiply_vector(&code.hsyn, &binary_error, &gf2_syn) == GF2_OK);
        for (size_t r = 0; r < 4; ++r) {
            unsigned parity = 0;
            for (size_t c = 0; c < 5; ++c) {
                Pauli p = stabilizer[r][c];
                unsigned x = p == PAULI_X || p == PAULI_Y;
                unsigned z = p == PAULI_Z || p == PAULI_Y;
                parity ^= (bits[c] & z) ^ (bits[c+5] & x);
            }
            CHECK(syn_data[r] == parity && gf2_data[r] == parity);
        }
    }
    quantum_code_free(&code);
    CHECK(code.hx.data == NULL && code.graph.edges == NULL);

    uint8_t ax[] = {1,0}, az[] = {0,1}; /* X and Z on same qubit. */
    BinaryMatrix bad_x = {2,1,ax}, bad_z = {2,1,az};
    CHECK(quantum_code_init(&code, &bad_x, &bad_z) == QUANTUM_NONCOMMUTING);
    CHECK(!code.hx.data && !code.graph.edges);
    bad_z.cols = 2;
    CHECK(quantum_code_init(&code, &bad_x, &bad_z) == QUANTUM_DIMENSION_MISMATCH);
    bad_z.cols = 1;
    ax[0] = 2;
    CHECK(quantum_code_init(&code, &bad_x, &bad_z) == QUANTUM_INVALID_ARGUMENT);
}

static void test_channel(void)
{
    PauliProbabilities prob;
    double prior[3];
    CHECK(quantum_depolarizing_probabilities(0.12, &prob) == QUANTUM_OK);
    CHECK(near(prob.p_i, 0.88) && near(prob.p_x, 0.04));
    CHECK(prob.p_x == prob.p_y && prob.p_y == prob.p_z);
    CHECK(quantum_prior_llr(&prob, prior) == QUANTUM_OK);
    for (size_t k = 0; k < 3; ++k) CHECK(near(prior[k], log(22)));
    prob = (PauliProbabilities){0.4,0.1,0.3,0.2};
    CHECK(quantum_prior_llr(&prob, prior) == QUANTUM_OK);
    CHECK(near(prior[0], log(4)) && near(prior[1], log(4.0/3)) && near(prior[2], log(2)));
    CHECK(quantum_depolarizing_probabilities(0, &prob) == QUANTUM_OK);
    CHECK(quantum_prior_llr(&prob, prior) == QUANTUM_INVALID_ARGUMENT);

    Pauli a[32], b[32], other[32];
    CHECK(quantum_depolarizing_sample(0.3, 2026, 7, 32, a) == QUANTUM_OK);
    CHECK(quantum_depolarizing_sample(0.3, 2026, 8, 32, other) == QUANTUM_OK);
    CHECK(memcmp(a, other, sizeof(a)) != 0);
    CHECK(quantum_depolarizing_sample(0.3, 2026, 7, 32, b) == QUANTUM_OK);
    CHECK(memcmp(a, b, sizeof(a)) == 0);
    CHECK(quantum_depolarizing_sample(0.3, 2027, 7, 32, other) == QUANTUM_OK);
    CHECK(memcmp(a, other, sizeof(a)) != 0);
    const char labels[] = "IXYZ";
    const char golden[] = "IIIZIZXIIIIIZIYIIXIIIIIIIIXIIIXX";
    printf("Depolarizing replay seed=2026 frame=7 p=0.3: ");
    for (size_t c = 0; c < 32; ++c) {
        CHECK(labels[a[c]] == golden[c]);
        putchar(labels[a[c]]);
    }
    putchar('\n');
    CHECK(quantum_depolarizing_sample(0, 0, 0, 32, a) == QUANTUM_OK);
    for (size_t c = 0; c < 32; ++c) CHECK(a[c] == PAULI_I);
    CHECK(quantum_depolarizing_sample(1, 0, 0, 32, a) == QUANTUM_OK);
    for (size_t c = 0; c < 32; ++c) CHECK(a[c] != PAULI_I);

    Pauli samples[16384];
    unsigned counts[4] = {0};
    CHECK(quantum_depolarizing_sample(0.3, 2026, 0, 16384, samples) == QUANTUM_OK);
    for (size_t c = 0; c < 16384; ++c) ++counts[samples[c]];
    CHECK(fabs(counts[0] - 0.7*16384) < 350);
    for (size_t k = 1; k < 4; ++k) CHECK(fabs(counts[k] - 0.1*16384) < 150);
    printf("Fixed sampler counts I/X/Y/Z: %u/%u/%u/%u\n", counts[0], counts[1], counts[2], counts[3]);
}

int main(void)
{
    test_pauli_and_graph();
    test_channel();
    printf("Quantum representation/channel passed (%u checks; 4^5 syndrome oracle)\n", checks);
    return EXIT_SUCCESS;
}
