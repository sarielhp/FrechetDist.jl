#!/usr/bin/env python3
"""Pairwise continuous Fréchet distances of the hydraulic pressure-sensor curves with
the Fred library (Rohde; pip package Fred-Frechet), for comparison with hydraulic.jl.

Curve k is sensor PSk: one vertex per machine cycle, its 6000 samples as coordinates
(Meintrup, Munteanu and Rohde, NeurIPS 2019). Only the first CYCLES vertices are used
(a comma-separated list of sizes, default 50,100; the full curves have 2205), and only
the first DIMS coordinates (default 6000, all of them). Writes
example_higher/output/hydraulic_fred_<CYCLES>x<DIMS>.csv per size, which hydraulic.jl
reads.

Fred finds the continuous Fréchet distance by binary search and stops once the bounds
are within `continuous_frechet_error` percent of the lower bound (default 1); it is
set to 1e-7 percent here, so hydraulic.jl can require agreement to 1e-8.

Python, not Ruby: Fred has Python bindings only. Usage, after
`ruby example_higher/download_hydraulic.rb` and installing Fred-Frechet (see README):
    python example_higher/fred_hydraulic.py [CYCLES[,CYCLES...] [DIMS]]
"""

import itertools
import os
import sys
import time

import numpy as np
import Fred

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data", "hydraulic")
OUTPUT = os.path.join(HERE, "output")
SENSORS = 6


def load_points(k, cycles, dims):
    path = os.path.join(DATA, f"PS{k}.txt")
    points = np.loadtxt(path, delimiter="\t", max_rows=cycles, dtype=np.float64)
    return points[:, :dims]


def run_size(points, cycles):
    curves = {k: Fred.Curve(np.ascontiguousarray(p[:cycles]), f"PS{k}")
              for k, p in points.items()}
    dim = curves[1].dimensions
    print(f"\nFred: {cycles} vertices per curve in R^{dim}")
    rows = []
    for i, j in itertools.combinations(range(1, SENSORS + 1), 2):
        start = time.perf_counter()
        result = Fred.continuous_frechet(curves[i], curves[j])
        seconds = time.perf_counter() - start
        rows.append((i, j, result.value, seconds))
        print(f"  PS{i}-PS{j}  {result.value:14.6f}  {seconds:9.3f} s  "
              f"({result.number_searches} decisions)")
    print(f"total {sum(r[3] for r in rows):.2f} s")
    os.makedirs(OUTPUT, exist_ok=True)
    path = os.path.join(OUTPUT, f"hydraulic_fred_{cycles}x{dim}.csv")
    with open(path, "w", encoding="utf-8") as io:
        io.write("i,j,distance,seconds\n")
        for i, j, value, seconds in rows:
            io.write(f"{i},{j},{value:.17g},{seconds:.6f}\n")
    print(f"written to {path}")


def main(args):
    sizes = [int(c) for c in args[0].split(",")] if args else [50, 100]
    dims = int(args[1]) if len(args) > 1 else 6000
    Fred.config.continuous_frechet_error = 1e-7
    Fred.config.number_threads = 1
    points = {k: load_points(k, max(sizes), dims) for k in range(1, SENSORS + 1)}
    for cycles in sizes:
        run_size(points, cycles)

if __name__ == "__main__":
    main(sys.argv[1:])
