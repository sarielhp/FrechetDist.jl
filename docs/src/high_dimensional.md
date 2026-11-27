# High-dimensional curves: comparison with Fred

This page compares FrechetDist with [Fred](https://pypi.org/project/Fred-Frechet/)
(Dennis Rohde's C++ library with Python bindings, version 1.14.5) on real curves in
ℝ^963, ℝ^6000 and ℝ^19200: the continuous Fréchet distance between every pair of
curves in three datasets. The scripts are in `example_higher/` of the repository; see
[Reproducing the results](@ref).

**In short:**

- **Faster.** On one thread, FrechetDist's exact algorithm takes 6.8 s for the 91
  traffic pairs (Fred: 110 s), 41 s for the 66 video pairs (Fred: 688 s), and 24
  minutes for the 15 hydraulic pairs (Fred: 6.3 hours): 16 to 17 times faster than
  Fred. On all 16 hardware threads of the machine it takes 1.3 s, 4.8 s and 160 s:
  **85, 143 and 142 times faster than Fred**.
- **Correct where Fred is not.** Fred returns a distance that is too low for 22 of
  the 172 pairs, by up to 3%, although it was asked for 10⁻⁹ relative accuracy. No
  matching of the two curves achieves the distance it returns. FrechetDist's results
  agree with an independent decision procedure.
- **Exact.** `frechet_c_compute_exact` has no tolerance: it returns the exact distance
  (a critical value of Alt and Godau, computed in floating point). `frechet_c_compute`
  returns the leash of an actual matching, with a certified lower bound, to a chosen
  tolerance; at the default 10⁻⁵ it takes 7.0 s and 24.7 s for the traffic and video
  sets, 9 and 14 times faster than Fred at the same tolerance.

## The datasets

| Set | Curves | Pairs | Vertices per curve | Dimension | Source |
|---|---|---|---|---|---|
| Freeway traffic | 14 days | 91 | 144 | 963 | [UCI PEMS-SF](https://archive.ics.uci.edu/dataset/204/pems+sf) |
| Video | 12 clips | 66 | 64–150 | 19 200 | [KTH human actions](https://www.csc.kth.se/cvap/actions/) |
| Hydraulic sensors | 6 sensors | 15 | 2205 | 6000 | [UCI hydraulic test rig](https://archive.ics.uci.edu/dataset/447/condition+monitoring+of+hydraulic+systems) |

- **Freeway traffic (PEMS-SF, CC BY 4.0).** Occupancy rates of 963 freeway sensors in
  the San Francisco Bay Area, every 10 minutes. Each day (2008–2009) is a curve: 144
  vertices, one per 10 minutes, each holding the 963 sensor readings. The first 14 days
  are used.
- **Video (KTH; Schuldt, Laptev and Caputo, ICPR 2004).** Grayscale 160×120 video at
  25 frames per second. Each frame is a vertex, its 19 200 pixels scaled to [0, 1]. A
  curve is the first repetition of an action in a video (64 to 150 frames), outdoors,
  for persons walking and boxing; the first 12 curves are used.
- **Hydraulic sensors (UCI).** Six pressure sensors PS1–PS6 of a hydraulic test rig.
  Each vertex is one machine cycle, its 6000 samples (60 s at 100 Hz) the coordinates.
  These are the curves of Meintrup, Munteanu and Rohde, [*Random Projections and
  Sampling Algorithms for Clustering of High-Dimensional Polygonal
  Curves*](https://arxiv.org/abs/1907.06969) (NeurIPS 2019), who report "roughly
  three hours" for a sequential implementation of Alt and Godau's algorithm.

## What is compared

| | Computes | Accuracy | Threads |
|---|---|---|---|
| `frechet_c_compute_exact(P, Q)` | the exact distance | exact: no tolerance | the tables use all threads; otherwise one |
| `frechet_c_compute(P, Q, true, tol)` | a matching (`Morphing`) whose leash is the distance, and a certified lower bound | relative `tol` (10⁻⁹ or the default 10⁻⁵) | one |
| Fred `continuous_frechet` | a value | asked for `continuous_frechet_error` percent (10⁻⁷ or 10⁻³, i.e. 10⁻⁹ or 10⁻⁵ relative) | one |

**`frechet_c_compute_exact`** implements the randomized algorithm of Har-Peled and
Raichel, "The Fréchet distance revisited and extended" (Section 6), for two curves.
The distance is one of the critical values of Alt and Godau (vertex–vertex,
vertex–segment and monotonicity events), the smallest one that the free-space decision
procedure accepts:

1. It tabulates, for every vertex of one curve and every segment of the other, the
   closest point on the segment's line and the squared distance to it, and all
   vertex–vertex distances: O(nmd) time, once, on all threads.
2. It draws 4N² random critical values from the tables (N the larger number of
   vertices), sorts them, and binary searches them with the decision procedure. This
   leaves an interval holding the distance and, with high probability, only
   O(N log N) critical values.
3. It lists the critical values in that interval. The monotonicity values come from a
   kinetic sweep along each segment, which keeps the intersection points of the
   growing spheres in a balanced search tree (O(log N) per event).
4. It binary searches the list. Each decision reads the tables, O(nm) time.

The result does not depend on the random choices, bit for bit; only the running time
does. `frechet_c_decide_exact(P, Q, r)`, the decision procedure (is the distance at
most r?), is exported too.

**`frechet_c_compute`** computes lower and upper bounds that meet: it simplifies the
curves, computes the VE-retractable Fréchet distance between the simplifications (a
lower bound) and a monotone matching (an upper bound), and refines until the two
agree to the tolerance. Its result is the leash of an actual matching, which it
returns, so it is never below the true distance.

**Fred** computes the continuous Fréchet distance by a binary search with the
Alt–Godau decision procedure, stopping when its bounds are within
`continuous_frechet_error` percent.

## Results

AMD Ryzen 7 5700G (8 cores, 16 hardware threads), 62 GB, Julia 1.13.3, Python 3.14.8,
Fred 1.14.5. Times are for all pairs of a set, excluding compilation (about 3 s, once
for all dimensions). "Wrong" counts the pairs more than 10⁻⁸ relative away from the
exact distance.

### High accuracy: 10⁻⁹

| Set | Fred | `frechet_c_compute` | exact, 1 thread | exact, 16 threads | Speed-up over Fred (exact, 1 / 16 threads) |
|---|---:|---:|---:|---:|---:|
| Traffic, ℝ^963, 91 pairs | 110.2 s | 26.9 s | 6.8 s | 1.3 s | 16× / 85× |
| Video, ℝ^19200, 66 pairs | 688.0 s | 54.1 s | 41.4 s | 4.8 s | 17× / 143× |
| Hydraulic, ℝ^6000, 15 pairs | 22 706 s (6.3 h) | 228.6 s | 1416 s | 159.9 s | 16× / 142× |

| Set | Fred: wrong | Fred: worst | `frechet_c_compute`: worst |
|---|---:|---:|---:|
| Traffic | 7 of 91 | −3.0% | +9.9·10⁻¹⁰ |
| Video | 13 of 66 | −1.5% | +9.9·10⁻¹⁰ |
| Hydraulic | 2 of 15 | −0.21% | −9.5·10⁻¹⁴ |

### Practical accuracy: 10⁻⁵

The package's default tolerance for `frechet_c_compute`, and Fred at 10⁻³ percent.
(Fred was not run at this tolerance on the hydraulic curves: about 3 hours.)

| Set | Fred | `frechet_c_compute` | exact, 1 thread | exact, 16 threads | Speed-up over Fred (`frechet_c_compute` / exact, 16 threads) |
|---|---:|---:|---:|---:|---:|
| Traffic, ℝ^963, 91 pairs | 60.4 s | 7.0 s | 6.8 s | 1.3 s | 8.6× / 46× |
| Video, ℝ^19200, 66 pairs | 350.1 s | 24.7 s | 41.4 s | 4.8 s | 14× / 73× |
| Hydraulic, ℝ^6000, 15 pairs | – | 227.5 s | 1416 s | 159.9 s | – |

- **`frechet_c_compute` keeps its tolerance:** every result is within 10⁻⁵ of the
  exact distance (worst +8.7·10⁻⁶), and never below it.
- **Fred does not:** the same 7 traffic and 13 video pairs as at 10⁻⁹ are too low by
  more than 10⁻⁵, by up to 3%.
- **The hydraulic curves** take `frechet_c_compute` the same time at both
  tolerances (227.5 s and 228.6 s), and at both its results agree with the exact
  distances to 10⁻¹³.

### Notes

- **Fred's errors are all too low.** The distance Fred returns for these pairs is not
  the leash of any matching: the decision procedure rejects it. For the hydraulic
  pair PS1–PS2, Fred returns 6083.53; the distance is 6096.29.
- **Two independent checks.** On the traffic and video sets the exact distances agree
  with `example_higher/check_frechet.jl`, a textbook Alt–Godau decision procedure
  written from scratch (no shared code), to 5·10⁻¹⁴ relative. On the hydraulic set,
  `frechet_c_compute` and `frechet_c_compute_exact` agree to 10⁻¹³.
- **Per pair.** On the hydraulic curves the exact algorithm takes 55–114 s per pair on
  one thread and 5.5–14.6 s on 16 threads, where Fred takes 18–37 minutes. With
  2205 vertices per curve, tabulating the vertex–segment distances (O(nmd)) is most
  of the work; it is the part that runs on all threads.
- **The paper's experiment.** Meintrup, Munteanu and Rohde report "roughly three
  hours" for their sequential Alt–Godau implementation on the hydraulic curves, and
  about 30 s with random projections and CUDA on Nvidia K40 GPUs. Their hardware and
  implementations differ, and they do not say whether the times are per distance or
  for all 15.

## Code

High-dimensional curves should use `Vector` points
(`GPntSeq{AnyDim,Float64,Vector{Float64}}`): the code is compiled once for every
dimension (static points, `Point{D}`, are compiled for each D, which takes minutes at
D = 6000). A matrix with one column per vertex converts directly:

```julia
using FrechetDist
using BasicCompGeometry: AnyDim, GPntSeq

const Curve = GPntSeq{AnyDim,Float64,Vector{Float64}}
curve(M::Matrix{Float64}) = Curve([M[:, j] for j in axes(M, 2)])  # D × n: one column per vertex

P, Q = curve(A), curve(B)

d = frechet_c_compute_exact(P, Q)        # the exact distance
frechet_c_decide_exact(P, Q, d)          # true: some matching has leash d

m = frechet_c_compute(P, Q)              # tolerance 1e-5 (the default)
m.leash, m.lower_bound                   # the distance is between these two
m = frechet_c_compute(P, Q, true, 1e-9)  # tighter
```

Start Julia with `-t auto` (or set `JULIA_NUM_THREADS=auto`) for the exact algorithm
to use all cores; the result is the same for any number of threads.

The same distance with Fred, as `example_higher/fred_curves.py` runs it:

```python
import Fred
Fred.config.continuous_frechet_error = 1e-7   # percent: 1e-9 relative
Fred.config.number_threads = 1
Fred.continuous_frechet(Fred.Curve(A), Fred.Curve(B)).value   # A, B: n × D arrays
```

## Reproducing the results

From the repository root (the examples have their own environment,
`example_higher/Project.toml`; see `example_higher/README.md` for installing Fred):

```sh
ruby example_higher/prepare_pems.rb                    # traffic: download, convert
ruby example_higher/prepare_kth.rb                     # video: download, convert (needs ffmpeg)
ruby example_higher/download_hydraulic.rb              # hydraulic sensors

python example_higher/fred_curves.py pems              # Fred: all pairs of days 1-14
python example_higher/fred_curves.py kth 1 12
python example_higher/fred_hydraulic.py 2205           # 6.3 hours

julia -t 1 example_higher/curves.jl pems               # frechet_c_compute, tolerance 1e-9
julia -t 1 example_higher/curves.jl kth 1 12
julia -t 1 example_higher/hydraulic.jl 2205

julia -t auto example_higher/exact_curves.jl pems      # frechet_c_compute_exact, and the
julia -t auto example_higher/exact_curves.jl kth 1 12  # errors of the others against it
julia -t auto example_higher/exact_curves.jl hydraulic

julia example_higher/check_frechet.jl pems 1 14        # independent check (slow)
```

The environment variable `FRECHET_TOL` (default `1e-9`) sets the relative tolerance of
`curves.jl`, `hydraulic.jl` and `fred_curves.py` (Fred: 100 × `FRECHET_TOL` percent);
results at other tolerances go to their own files.
