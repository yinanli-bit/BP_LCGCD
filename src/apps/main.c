#include "bplcgcd/decoder.h"
#include "bplcgcd/simulation.h"
#include "bplcgcd/slva.h"

#include <ctype.h>
#include <errno.h>
#include <inttypes.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

enum { CLI_OK = 0, CLI_USAGE_ERROR = 2, CLI_NOT_IMPLEMENTED = 3,
       CLI_FAILURE = 4, CLI_TRUNCATED = 5 };

static void print_help(void)
{
    puts("BP_LCGCD - single-frame GCD, bounded simulation and SLVA demo\n"
         "Usage: bp_lcgcd --help\n"
         "       bp_lcgcd --demo-slva\n"
         "       bp_lcgcd --decoder gcd --matrix H.txt [--mode MODE] INPUT [OPTIONS]\n"
         "       bp_lcgcd --decoder gcd --matrix H.txt --simulate --ebn0 DB\n"
         "         --frames N --errors N --csv PATH [--seed N] [--limit N]\n"
         "\n"
         "Modes and INPUT (default: classical):\n"
         "  classical     --llr \"L0 L1 ...\"\n"
         "  syndrome      --syndrome BITS --weights \"w0 w1 ...\" [--hint BITS]\n"
         "  bp-syndrome   --syndrome BITS --llr \"L0 L1 ...\" (BP posterior input only)\n"
         "\n"
         "Options:\n"
         "  -h, --help       Show this help (also shown with no arguments).\n"
         "  --demo-slva      Run the fixed S10-S13 example (stage-order local paths).\n"
         "  --decoder NAME   gcd, lcgcd, bp, nms, bp-lcgcd, bp-osd; only gcd is implemented.\n"
         "  --matrix PATH    H: positive rows/cols followed by whitespace-separated bits.\n"
         "  --mode MODE      Select one of the modes above.\n"
         "  --llr VALUES     Exactly N finite numbers, quoted and whitespace-separated.\n"
         "  --weights VALUES Exactly N finite nonnegative flip costs, quoted.\n"
         "  --syndrome BITS  Contiguous bits in original H row order.\n"
         "  --hint BITS      Contiguous reference bits in original column order (default zero).\n"
         "  --limit N        Complete candidate budget; 0=unlimited (default).\n"
         "  --trace          Print preparation, every evaluated candidate and stop condition.\n"
         "  --simulate       All-zero classical BPSK/AWGN at one Eb/N0 point.\n"
         "  --ebn0 DB        Finite Eb/N0 in dB for simulation.\n"
         "  --frames N       Positive simulation frame cap.\n"
         "  --errors N       Positive decoded frame-error cap.\n"
         "  --seed N         Simulation seed (default 12344567).\n"
         "  --csv PATH       Write CSV header + summary; codeword BER uses N bits.\n"
         "\n"
         "Decoder input/output bits use original order, index 0 first. Trace f_I/f_P use\n"
         "systematic order and print the zero-based sys_to_orig mapping.\n"
         "Demo paths use stage order (f_I[A],f_P).\n"
         "Exit codes: 0=help/proven optimum, 2=usage error, 3=not implemented,\n"
         "            4=input/decode/I/O failure, 5=one or more truncated frames.");
}

static int parse_soft(const char *text, double *values, size_t length, bool nonnegative)
{
    for (size_t j = 0; j < length; ++j) {
        char *end;
        errno = 0;
        values[j] = strtod(text, &end);
        if (end == text || errno == ERANGE || !isfinite(values[j]) ||
            (nonnegative && values[j] < 0) || (*end && !isspace((unsigned char)*end))) return 0;
        text = end;
    }
    while (isspace((unsigned char)*text)) ++text;
    return *text == '\0';
}

static int parse_bits(const char *text, BitVector *bits)
{
    if (strlen(text) != bits->length) return 0;
    for (size_t j = 0; j < bits->length; ++j) {
        if (text[j] != '0' && text[j] != '1') return 0;
        bits->data[j] = (uint8_t)(text[j] - '0');
    }
    return 1;
}

static const char *status_name(DecodeStatus status)
{
    static const char *const names[] = {"NOT_RUN", "SUCCESS", "TRUNCATED", "NOT_IMPLEMENTED",
        "INVALID_ARGUMENT", "DIMENSION_MISMATCH", "INVALID_CANDIDATE", "NUMERIC_ERROR",
        "ALLOCATION_FAILED", "INCONSISTENT_SYNDROME"};
    return names[status];
}

static int parse_count(const char *text, uint64_t *value)
{
    char *end;
    errno = 0;
    uintmax_t parsed = strtoumax(text, &end, 10);
    if (*text < '0' || *text > '9' || *end || errno == ERANGE || parsed > UINT64_MAX) return 0;
    *value = (uint64_t)parsed;
    return 1;
}

/* Fixed teaching example: edit inputs here; S12 runs inside slva_init. */
static int run_slva_demo(void)
{
    enum { ROWS=3, N=6, DELTA=2, K=4, L=DELTA+K };
    uint8_t h_data[ROWS*N] = {
        0,1,1,0,1,0,
        0,0,1,1,0,1,
        0,1,0,1,1,1
    };
    double flip_weights[N] = {1000,3,6,9,8,7}; /* Original column order. */
    uint8_t residual_data[ROWS] = {1,1,0}; /* Pure syndrome, zero reference. */
    BinaryMatrix h = {ROWS,N,h_data};
    BitVector s_res = {ROWS,residual_data};
    BinarySystematic form = {0};
    SlvaGenerator list = {0};
    size_t basis_order[DELTA], masks[L], columns[L];
    uint8_t pa[DELTA*K], sa[DELTA], scratch_data[ROWS], path_data[L];
    LocalConstraints local = {.basis_order=basis_order,
        .p_a={DELTA,K,pa}, .p_ac={0,K,NULL}, .s_a={DELTA,sa}, .s_ac={0,NULL}};
    BitVector scratch = {ROWS,scratch_data}, path = {L,path_data};
    SyndromeTrellis trellis = {.length=L,.column_masks=masks,.original_columns=columns};
    SlvaLimits limits = {.max_candidates=128,.max_bytes=64*1024,.path_limit=0};
    int exit_code = CLI_FAILURE;
    const char *stage = "S03 systematize";
    if (binary_matrix_systematize(&h,&form) != GF2_OK) goto cleanup;
    stage = "S10 local_constraints_prepare";
    if (local_constraints_prepare(&form,flip_weights,N,&s_res,DELTA,&local,&scratch) != GF2_OK) goto cleanup;
    stage = "S11 syndrome_trellis_prepare";
    if (syndrome_trellis_prepare(&form,&local,flip_weights,N,4,&trellis) != TRELLIS_OK) goto cleanup;
    printf("S10-S13 demo: N=%d rows=%d m=%zu K=%d delta=%d\n",N,ROWS,form.rank,K,DELTA);
    printf("trellis length=%zu states=%zu start=0 target=%zu\n",
           trellis.length,trellis.state_count,trellis.target);
    stage = "S13 slva_init / S12 viterbi_decode";
    if (slva_init(&list,&trellis,&limits) != SLVA_OK) goto cleanup;
    stage = "S13 slva_next";
    SlvaStatus status;
    double cost;
    while ((status=slva_next(&list,&path,&cost)) == SLVA_OK) {
        printf("rank=%" PRIu64 " path(stage)=",list.paths_returned);
        for (size_t j = 0; j < path.length; ++j) putchar('0'+path.data[j]);
        printf(" cost=%.17g\n",cost);
    }
    if (status != SLVA_DONE) goto cleanup;
    printf("SLVA demo DONE paths=%" PRIu64 "\n",list.paths_returned);
    exit_code = CLI_OK;
cleanup:
    if (exit_code != CLI_OK) fprintf(stderr,"DEMO_ERROR: %s failed.\n",stage);
    slva_free(&list);
    binary_systematic_free(&form);
    return exit_code;
}

int main(int argc, char *argv[])
{
    if (argc == 1 || (argc == 2 && (!strcmp(argv[1], "--help") || !strcmp(argv[1], "-h")))) {
        print_help();
        return CLI_OK;
    }
    if (argc == 2 && !strcmp(argv[1], "--demo-slva")) return run_slva_demo();
    const char *decoder = NULL, *matrix_path = NULL, *mode = NULL, *llr_text = NULL;
    const char *weights_text = NULL, *syndrome_text = NULL, *hint_text = NULL, *limit_text = NULL;
    const char *ebn0_text = NULL, *frames_text = NULL, *errors_text = NULL, *seed_text = NULL, *csv_path = NULL;
    const struct { const char *name; const char **value; } options[] = {
        {"--decoder", &decoder}, {"--matrix", &matrix_path}, {"--mode", &mode},
        {"--llr", &llr_text}, {"--weights", &weights_text}, {"--syndrome", &syndrome_text},
        {"--hint", &hint_text}, {"--limit", &limit_text}, {"--ebn0", &ebn0_text},
        {"--frames", &frames_text}, {"--errors", &errors_text}, {"--seed", &seed_text}, {"--csv", &csv_path}
    };
    bool trace = false, simulate = false;
    for (int i = 1; i < argc; ++i) {
        if (!strcmp(argv[i], "--trace") && !trace) { trace = true; continue; }
        if (!strcmp(argv[i], "--simulate") && !simulate) { simulate = true; continue; }
        const char **slot = NULL;
        for (size_t j = 0; j < sizeof(options) / sizeof(*options); ++j) {
            if (!strcmp(argv[i], options[j].name)) { slot = options[j].value; break; }
        }
        if (!slot || *slot || i + 1 == argc) {
            fputs("USAGE_ERROR: unknown/duplicate option or missing value. Use --help.\n", stderr);
            return CLI_USAGE_ERROR;
        }
        *slot = argv[++i];
    }
    if (!decoder || !decoder_is_known(decoder)) {
        fputs("USAGE_ERROR: unknown decoder or missing --decoder. Use --help.\n", stderr);
        return CLI_USAGE_ERROR;
    }
    if (strcmp(decoder, "gcd")) {
        fprintf(stderr, "NOT_IMPLEMENTED: decoder '%s' is not available.\n", decoder);
        return CLI_NOT_IMPLEMENTED;
    }
    if (!mode) mode = "classical";
    DecodeInput input = {0};
    if (!strcmp(mode, "classical")) input.mode = DECODE_CLASSICAL;
    else if (!strcmp(mode, "syndrome")) input.mode = DECODE_SYNDROME;
    else if (!strcmp(mode, "bp-syndrome")) input.mode = DECODE_BP_SYNDROME;
    else { fputs("USAGE_ERROR: unknown mode. Use --help.\n", stderr); return CLI_USAGE_ERROR; }
    bool valid_options;
    if (simulate) valid_options = input.mode == DECODE_CLASSICAL && !llr_text && !weights_text &&
                                 !syndrome_text && !hint_text && ebn0_text && frames_text && errors_text && csv_path;
    else if (input.mode == DECODE_CLASSICAL) valid_options = llr_text && !weights_text && !syndrome_text && !hint_text;
    else if (input.mode == DECODE_SYNDROME) valid_options = syndrome_text && weights_text && !llr_text;
    else valid_options = syndrome_text && llr_text && !weights_text && !hint_text;
    if (!matrix_path || !valid_options || (!simulate && (ebn0_text || frames_text || errors_text || seed_text || csv_path))) {
        fputs("USAGE_ERROR: missing matrix/input or incompatible mode options. Use --help.\n", stderr);
        return CLI_USAGE_ERROR;
    }
    DecodeConfig config = {.decoder_name=decoder, .trace=trace ? stdout : NULL};
    if (limit_text && !parse_count(limit_text, &config.candidate_limit)) {
        fputs("USAGE_ERROR: --limit requires a nonnegative integer.\n", stderr);
        return CLI_USAGE_ERROR;
    }
    SimulationConfig sim = {.seed=12344567};
    if (simulate && (!parse_soft(ebn0_text, &sim.ebn0_db, 1, false) ||
        !parse_count(frames_text, &sim.max_frames) || !sim.max_frames ||
        !parse_count(errors_text, &sim.max_frame_errors) || !sim.max_frame_errors ||
        (seed_text && !parse_count(seed_text, &sim.seed)))) {
        fputs("USAGE_ERROR: simulation requires finite Eb/N0, positive caps and a nonnegative seed.\n", stderr);
        return CLI_USAGE_ERROR;
    }

    BinaryMatrix h = {0};
    BinarySystematic form = {0};
    BitVector target = {0}, hint = {0};
    DecodeResult result = {0};
    double *soft = NULL;
    int exit_code = CLI_FAILURE;
    FILE *stream = fopen(matrix_path, "r");
    if (!stream) { fprintf(stderr, "INPUT_ERROR: cannot open matrix '%s'.\n", matrix_path); goto cleanup; }
    Gf2Status loaded = binary_matrix_load(stream, &h);
    fclose(stream);
    if (loaded != GF2_OK) { fputs("INPUT_ERROR: cannot load matrix (format, I/O or allocation).\n", stderr); goto cleanup; }
    if (simulate) {
        if (binary_matrix_systematize(&h, &form) != GF2_OK) {
            fputs("INPUT_ERROR: offline systematization failed.\n", stderr);
            goto cleanup;
        }
        SimulationResult summary = {0};
        DecodeStatus status = simulation_run(&h, &form, &config, &sim, &summary);
        if (status != DECODE_SUCCESS) {
            fprintf(stderr, "DECODE_ERROR: simulation %s\n", status_name(status));
            goto cleanup;
        }
        FILE *csv = fopen(csv_path, "w");
        if (!csv) { fprintf(stderr, "INPUT_ERROR: cannot open CSV '%s'.\n", csv_path); goto cleanup; }
        int written = simulation_write_csv(csv, &config, &sim, &summary);
        if (fclose(csv) != 0) written = 0;
        if (!written) { fputs("IO_ERROR: cannot write CSV summary.\n", stderr); goto cleanup; }
        printf("simulation frames=%" PRIu64 " frame_errors=%" PRIu64 " truncated_frames=%" PRIu64
               " stop_reason=%s csv=%s\n", summary.frames, summary.frame_errors, summary.truncated_frames,
               summary.error_limit_hit ? "error_limit" : "frame_limit", csv_path);
        exit_code = summary.truncated_frames ? CLI_TRUNCATED : CLI_OK;
        goto cleanup;
    }
    if (h.cols > SIZE_MAX / sizeof(*soft)) goto allocation_error;
    soft = malloc(h.cols * sizeof(*soft));
    if (!soft || bit_vector_init(&result.estimate, h.cols) != GF2_OK ||
        (syndrome_text && bit_vector_init(&target, h.rows) != GF2_OK) ||
        (hint_text && bit_vector_init(&hint, h.cols) != GF2_OK)) goto allocation_error;
    if (!parse_soft(llr_text ? llr_text : weights_text, soft, h.cols, input.mode == DECODE_SYNDROME) ||
        (syndrome_text && !parse_bits(syndrome_text, &target)) ||
        (hint_text && !parse_bits(hint_text, &hint))) {
        fputs("USAGE_ERROR: invalid values/bits or length mismatch with H.\n", stderr);
        exit_code = CLI_USAGE_ERROR;
        goto cleanup;
    }
    if (binary_matrix_systematize(&h, &form) != GF2_OK) {
        fputs("INPUT_ERROR: offline systematization failed.\n", stderr);
        goto cleanup;
    }
    input.h = &h;
    input.length = h.cols;
    input.systematic = &form;
    if (input.mode == DECODE_SYNDROME) input.weights = soft;
    else input.llr = soft;
    if (syndrome_text) input.s_target = &target;
    if (hint_text) input.z_hint = &hint;
    DecodeStatus status = decoder_decode(&config, &input, &result);
    if (status != DECODE_SUCCESS && status != DECODE_TRUNCATED) {
        fprintf(stderr, "DECODE_ERROR: %s\n", status_name(status));
        goto cleanup;
    }
    static const char *const stops[] = {"none", "bound", "exhausted", "budget"};
    printf("mode=%s status=%s candidate_valid=%d optimality_proven=%d stop_reason=%s\n",
           mode, status_name(status), (int)result.candidate_valid, (int)result.optimality_proven,
           stops[result.stop_reason]);
    if (result.candidate_valid) {
        printf("estimate_type=%s estimate=", input.mode == DECODE_CLASSICAL ? "codeword" : "error");
        for (size_t j = 0; j < h.cols; ++j) putchar('0' + result.estimate.data[j]);
        printf(" cost=%.17g\n", result.cost);
    }
    printf("candidates_popped=%" PRIu64 " candidates_evaluated=%" PRIu64
           " nodes_expanded=%" PRIu64 " batches=%" PRIu64 "\n", result.stats.candidates_popped,
           result.stats.candidates_evaluated, result.stats.nodes_expanded, result.stats.batches);
    exit_code = status == DECODE_TRUNCATED ? CLI_TRUNCATED : CLI_OK;
    goto cleanup;
allocation_error:
    fputs("INPUT_ERROR: allocation failed or dimensions too large.\n", stderr);
cleanup:
    free(soft);
    bit_vector_free(&result.estimate);
    bit_vector_free(&hint);
    bit_vector_free(&target);
    binary_systematic_free(&form);
    binary_matrix_free(&h);
    return exit_code;
}
