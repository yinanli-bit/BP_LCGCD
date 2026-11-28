"""Check the direct main.c simulation against TASK2's deterministic results."""
import csv
from pathlib import Path
import subprocess
import sys

exe, root, output = map(Path, sys.argv[1:4])
output.mkdir(parents=True, exist_ok=True)
hx = root / "assets/quantum/five_qubit_hx.txt"
hz = root / "assets/quantum/five_qubit_hz.txt"
fields = ["p", "decoder", "frames", "logical_failures", "ler", "bp_iterations_total",
          "bp_iterations_avg", "lc_searches_total", "lc_searches_avg", "bp_converged_frames"]

for option in ([], ["--help"], ["-h"]):
    help_result = subprocess.run([str(exe)] + option, capture_output=True, text=True)
    assert help_result.returncode == 0 and "--quantum" in help_result.stdout
    assert not any(old in help_result.stdout for old in ("--demo-slva", "--decoder", "--simulate", "--compare"))
removed = subprocess.run([str(exe), "--demo-slva"], capture_output=True)
assert removed.returncode == 2


def run(name, start, stop, step, frames=64, failures=0, iterations=1, limit=0):
    path = output / f"{name}.csv"
    command = [str(exe), "--quantum", str(hx), str(hz), str(start), str(stop), str(step),
               str(frames), str(failures), str(iterations), ".625", "2", str(limit), "2026", str(path)]
    result = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True, timeout=20)
    assert result.returncode == 0, result.stderr
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        assert reader.fieldnames == fields
        rows = list(reader)
    for row in rows:
        n = int(row["frames"])
        assert float(row["ler"]) == int(row["logical_failures"]) / n
        assert float(row["bp_iterations_avg"]) == int(row["bp_iterations_total"]) / n
        assert float(row["lc_searches_avg"]) == int(row["lc_searches_total"]) / n
    return rows, path


golden, path = run("golden", .12, .12, .04)
assert [row["decoder"] for row in golden] == ["QNMS", "BP-LCGCD"]
assert [int(row["logical_failures"]) for row in golden] == [23, 10]
assert [int(row["bp_iterations_total"]) for row in golden] == [64, 64]
assert [int(row["bp_converged_frames"]) for row in golden] == [43, 43]
assert [int(row["lc_searches_total"]) for row in golden] == [0, 172]
assert int(golden[0]["logical_failures"]) > 64 - 43  # Syndrome match is not logical success.
_, replay = run("replay", .12, .12, .04)
assert path.read_bytes() == replay.read_bytes()

multi, _ = run("multi", .04, .12, .04)
assert len(multi) == 6
for point in range(3):
    pair = multi[2*point:2*point+2]
    assert abs(float(pair[0]["p"]) - (.04 + point*.04)) < 1e-15
    assert pair[0]["frames"] == pair[1]["frames"] == "64"
    single, _ = run(f"single{point}", pair[0]["p"], pair[0]["p"], .04)
    assert pair == single

stopped, _ = run("stopped", .12, .12, .04, failures=2)
assert stopped[0]["frames"] == stopped[1]["frames"]
assert int(stopped[0]["frames"]) < 64 and all(int(row["logical_failures"]) >= 2 for row in stopped)
capped, _ = run("capped", .12, .12, .04, limit=1)
assert 0 < int(capped[1]["lc_searches_total"]) <= int(capped[1]["frames"])

bad = subprocess.run([str(exe), "--quantum", str(hx), str(hz), ".1", ".2", "0",
                      "64", "0", "1", ".625", "2", "0", "2026", str(output / "bad.csv")],
                     stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
assert bad.returncode == 2
print("Quantum CLI passed: TASK2 reference, shared frame counts, point loop, replay, failure stop and candidate cap")
