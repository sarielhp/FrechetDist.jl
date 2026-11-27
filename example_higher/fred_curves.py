#!/usr/bin/env python3
"""Pairwise continuous Fréchet distances among curves of a prepared set with the Fred
library (Rohde; pip package Fred-Frechet), for comparison with curves.jl.

Curves FIRST to LAST of data/curves/SET/index.csv (written by prepare_pems.rb or
prepare_kth.rb; default 1 to 14). Each curve file holds dims x vertices Float64
values, vertex after vertex. Writes example_higher/output/SET_FIRST-LAST_fred.csv,
which curves.jl reads. Fred's binary search stops within `continuous_frechet_error`
percent: 100 times the relative tolerance FRECHET_TOL (environment, default 1e-9, as
for curves.jl; results at other tolerances go to files with _tolFRECHET_TOL in their
names). One thread.

Python, not Ruby: Fred has Python bindings only. Usage:
    python example_higher/fred_curves.py SET [FIRST [LAST]]
"""

import csv
import itertools
import os
import sys
import time

import numpy as np
import Fred

HERE = os.path.dirname(os.path.abspath(__file__))
CURVES = os.path.join(HERE, "data", "curves")
OUTPUT = os.path.join(HERE, "output")


def read_index(name):
    with open(os.path.join(CURVES, name, "index.csv"), encoding="utf-8") as io:
        return list(csv.DictReader(io))


def read_curve(name, row):
    dims, vertices = int(row["dims"]), int(row["vertices"])
    values = np.fromfile(os.path.join(CURVES, name, row["file"]), dtype="<f8")
    # Column-major dims x vertices is row-major vertices x dims: one row per vertex.
    points = np.ascontiguousarray(values.reshape(vertices, dims))
    return Fred.Curve(points, row["name"])


def main(args):
    if not args:
        sys.exit("usage: fred_curves.py SET [FIRST [LAST]]")
    name = args[0]
    index = read_index(name)
    first = int(args[1]) if len(args) > 1 else 1
    last = int(args[2]) if len(args) > 2 else min(14, len(index))
    rows = index[first - 1:last]
    tol_env = os.environ.get("FRECHET_TOL", "1e-9")
    tol = float(tol_env)
    suffix = "" if tol == 1e-9 else f"_tol{tol_env}"
    # Percent, rounded to the decimal value it stands for (1e-7, not 1.0000000000000001e-07).
    Fred.config.continuous_frechet_error = float(f"{100 * tol:.15g}")
    Fred.config.number_threads = 1
    curves = [read_curve(name, r) for r in rows]
    print(f"Fred: {name} curves {first}-{last} in R^{curves[0].dimensions}")
    results = []
    for i, j in itertools.combinations(range(len(curves)), 2):
        start = time.perf_counter()
        value = Fred.continuous_frechet(curves[i], curves[j]).value
        seconds = time.perf_counter() - start
        results.append((rows[i]["name"], rows[j]["name"], value, seconds))
    print(f"{len(results)} pairs: {sum(r[3] for r in results):.2f} s")
    os.makedirs(OUTPUT, exist_ok=True)
    path = os.path.join(OUTPUT, f"{name}_{first}-{last}{suffix}_fred.csv")
    with open(path, "w", encoding="utf-8") as io:
        io.write("a,b,distance,seconds\n")
        for a, b, value, seconds in results:
            io.write(f"{a},{b},{value:.17g},{seconds:.6f}\n")
    print(f"written to {path}")


if __name__ == "__main__":
    main(sys.argv[1:])
