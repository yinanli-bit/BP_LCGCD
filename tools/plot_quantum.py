#!/usr/bin/env python3
"""Plot the three quantum comparison curves from the simulation CSV."""
import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("csv", type=Path)
parser.add_argument("--output-dir", type=Path)
args = parser.parse_args()
with args.csv.open(newline="") as stream:
    rows = list(csv.DictReader(stream))
if not rows:
    parser.error("CSV has no simulation rows")
output = args.output_dir or args.csv.parent / "plots"
output.mkdir(parents=True, exist_ok=True)

for field, ylabel, filename in (
    ("ler", "Logical error rate", "ler.png"),
    ("bp_iterations_avg", "Average BP iterations per frame", "bp_iterations.png"),
    ("lc_searches_avg", "Average LC-GCD candidates per frame", "lc_searches.png"),
):
    fig, ax = plt.subplots()
    for decoder in dict.fromkeys(row["decoder"] for row in rows):
        points = sorted((row for row in rows if row["decoder"] == decoder), key=lambda row: float(row["p"]))
        ax.plot([float(row["p"]) for row in points], [float(row[field]) for row in points], "o-", label=decoder)
    ax.set(xlabel="Physical Pauli error rate p", ylabel=ylabel)
    ax.set_ylim(bottom=0)
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    path = output / filename
    fig.savefig(path)
    plt.close(fig)
    print(path)
