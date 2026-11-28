#!/usr/bin/env bash
set -e
cd -- "$(dirname -- "$0")"

# Edit these values, then run ./run_quantum_sim.sh.
EXECUTABLE=./build/quantum/bp_lcgcd
HX=assets/quantum/five_qubit_hx.txt
HZ=assets/quantum/five_qubit_hz.txt
P_START=0.04
P_STOP=0.12
P_STEP=0.04
MAX_FRAMES=1000
MAX_LOGICAL_FAILURES=100
BP_MAX_ITERATIONS=12
NMS_SCALE=0.625
DELTA=2
CANDIDATE_LIMIT=128
SEED=2026
OUTPUT=results/quantum/run.csv

if [[ ! -x "$EXECUTABLE" ]]; then
    echo 'Build first: cmake -S . -B build/quantum -DCMAKE_BUILD_TYPE=Release && cmake --build build/quantum -j 4' >&2
    exit 2
fi
mkdir -p -- "$(dirname -- "$OUTPUT")"
# TASK1/2 teaching prints stay unchanged; point summaries go to stderr.
exec "$EXECUTABLE" --quantum "$HX" "$HZ" "$P_START" "$P_STOP" "$P_STEP" \
    "$MAX_FRAMES" "$MAX_LOGICAL_FAILURES" "$BP_MAX_ITERATIONS" "$NMS_SCALE" \
    "$DELTA" "$CANDIDATE_LIMIT" "$SEED" "$OUTPUT" > /dev/null
