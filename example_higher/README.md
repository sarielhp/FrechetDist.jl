# Fréchet distance between curves in higher dimensions

These examples treat human motion as curves in ℝ^46 and compare them with the package.
The data comes from the [CMU Graphics Lab Motion Capture Database](http://mocap.cs.cmu.edu/).
Each trial is a skeleton recorded at 120 Hz. Every frame becomes one vertex of a polygonal
curve, and its coordinates are the skeleton's joint angles in degrees.

## Running

From the repository root:

```sh
ruby example_higher/download.rb              # 41 trials of 5 subjects, about 9 MB, into example_higher/data/
julia example_higher/compare_two.jl          # one walk vs. one run, under five distances
julia example_higher/compare_two.jl 09_01 16_35 4   # any two trials, every 4th frame
julia example_higher/distance_matrix.jl      # all 820 pairs (about 45 s), nearest neighbours
julia example_higher/distance_matrix.jl 4 1.05      # every 4th frame, 1.05-approximation
```

The scripts start with `using QuickEnv # local`, which activates this directory's
`Project.toml`; FrechetDist resolves to the working tree. Both `data/` and `output/`
are ignored by git.

## Files

| File | Contents |
|---|---|
| `trials.txt` | The trials used, each with an activity class (walk, run, jump) taken from the database's description |
| `download.rb` | Fetches each subject's skeleton (`.asf`) and each trial's motion (`.amc`); files already present are kept |
| `mocap.jl` | Reads `.amc` files into `PntSeq{46,Float64}` curves |
| `compare_two.jl` | Exact and approximate Fréchet, VE-retractable Fréchet, discrete Fréchet and DTW between two trials, and the order the first four must satisfy |
| `distance_matrix.jl` | All-pairs exact Fréchet distances, nearest neighbours, mean distance between activities, and `output/distances.csv` |
| `point_types.jl` | Static against vector points, per dimension: speed of the point kernels and compile time (see [below](#static-or-vector-points)) |

## The curves

- **The root bone is dropped.** Its values are the position in the lab and the heading, so they depend on where the subject stood and which way it faced, not on the motion.
- **The hands, fingers and thumbs are dropped.** The database tracks them poorly, and their angles flip by 150 to 180 degrees between frames.
- **Frames with missing data are dropped.** These are frames where a bone reads exactly `0 0 0`; at most 7 frames per trial are affected.

The remaining 46 angles change by at most about 47 degrees between consecutive frames.

The Fréchet distance takes a maximum, so a single bad frame can decide it. That is why the
data is cleaned at all.

## What to expect

With exact distances, the nearest other trial has the same activity for 40 of the 41
trials. Excluding trials of the same subject, the nearest trial matches for 35 of 41.
The misses are mostly jumps matched with walks: only two subjects contribute jumps, and
a jump trial is mostly standing and stepping.

Keeping every 4th frame with a 1.05-approximation (`distance_matrix.jl 4 1.05`) takes about
7 s instead of 45 s, and gives nearly the same results (40/41 and 36/41).

## Comparison with a published experiment: curves in ℝ⁶⁰⁰⁰

Meintrup, Munteanu and Rohde, [*Random Projections and Sampling Algorithms for Clustering
of High-Dimensional Polygonal Curves*](https://arxiv.org/abs/1907.06969) (NeurIPS 2019,
Section 5), compute continuous Fréchet distances between six curves taken from the
[UCI hydraulic test rig data](https://archive.ics.uci.edu/dataset/447/condition+monitoring+of+hydraulic+systems).
Curve k is pressure sensor PSk. It has one vertex per machine cycle, 2205 in all, and
each vertex holds that cycle's 6000 samples. They report "roughly three hours" for
their sequential Alt–Godau implementation. Random projection and a CUDA version
together bring this to about 30 s "on average", on a cluster with Nvidia K40 GPUs. The
paper does not say whether these times are per distance or for all 15, and its code
link (dennisrohde.work/rp4frechet-code) no longer responds, so the comparison is
approximate.

```sh
ruby example_higher/download_hydraulic.rb               # 77 MB zip, PS1..PS6 extracted (450 MB)
python example_higher/fred_hydraulic.py [CYCLES[,CYCLES...] [DIMS]]  # Fred (Rohde's library), optional
julia example_higher/hydraulic.jl [CYCLES[,CYCLES...] [DIMS]]        # this package; lists Fred's results if present
```

Both scripts use only the first CYCLES vertices of each curve (a list of sizes,
default `50,100`) and the first DIMS coordinates (default 6000). A short run therefore
checks agreement and measures how the time scales. `hydraulic.jl` compiles once for all
sizes. It estimates the time for the full curves from the growth between the last two
sizes, and fails if a distance differs from Fred's by more than 10⁻⁸ relative.
Running with `2205 6000` uses the paper's curves; the hardware and implementations
differ.

[Fred](https://pypi.org/project/Fred-Frechet/) needs a C++ compiler; it has no wheel
for recent Pythons. Install it into the project cache:

```sh
python3 -m venv ~/.cache/FrechetDist/fred-venv
~/.cache/FrechetDist/fred-venv/bin/pip install Fred-Frechet
~/.cache/FrechetDist/fred-venv/bin/python example_higher/fred_hydraulic.py 50
```

- **Fred's tolerance:** Fred computes the continuous Fréchet distance by binary search, to within `continuous_frechet_error` percent (default 1%). The script sets this to 10⁻⁷ percent and uses one thread.
- **Fred's scaling:** with all 6000 coordinates, Fred takes 1.7 s, 7.1 s and 29.3 s for all 15 pairs at 50, 100 and 200 vertices. That is quadratic growth, or about an hour for the full curves. 42 of these 45 distances needed no search: Fred's quick lower and upper bounds already met. These curves are easy instances.
- **Compile time:** the script uses `Vector` points (`GPntSeq{AnyDim}`), so the code is compiled once for every dimension, in about 3 s. With static points (`SVector{D}`) it was compiled for each dimension: about 25 s for D = 1000, 100 s for D = 2000, and 11–14 minutes for D = 6000. The times the script reports exclude it.

### Results so far (truncated curves, all 6000 coordinates, one thread)

| Vertices per curve | This package, `Vector` points | static points | Fred | Largest relative difference |
|---|---|---|---|---|
| 50  | 0.50 s | 3.45 s | 1.68 s  | 6.8·10⁻¹⁰ |
| 100 | 0.87 s | 5.55 s | 6.81 s  | 8.2·10⁻¹⁰ |
| 200 | 1.63 s | 8.81 s | 27.77 s | 8.6·10⁻¹⁰ |

Times are for all 15 pairs; the differences to Fred are those of the `Vector`-point run (the static run gave the same).

- **Scaling:** Fred grows quadratically. This package (`Vector` points) grew as about n^0.9 from 100 to 200 vertices, which would put the full curves near 15 s against about an hour for Fred. That extrapolation rests on two sizes and is not a measurement.
- **The distances do not change with the size:** all 15 are the same at 50, 100 and 200 vertices, so the early cycles decide them. The full run (`2205 6000`) is the real test.

## Static or vector points

Points are static vectors (`SVector{D}`, BasicCompGeometry's `Point{D}`). Their
length is part of the type, so the compiler unrolls the loops over the coordinates,
and points live on the stack or inline in arrays. That is ideal in the plane. In high
dimensions it backfires: the code is compiled again for every D, and both the compile
time and the unrolled code grow with D. A `Vector{Float64}` point has its length
known only at run time. Its code is compiled once for all dimensions and runs a
plain vectorized loop, but every point it creates is a heap allocation.
BasicCompGeometry's `GPntSeq{AnyDim,T,Vector{T}}` holds such points, and every
FrechetDist algorithm takes either kind (`src/fp_geometry_vec.jl` has the loop kernels).

`point_types.jl` measures where the crossover lies, for the three kernels the Fréchet
code spends its time in. `dist` is the distance between two points. `dps` is the
distance from a point to a segment, the main work of the free-space computations.
`comb` creates a point on a segment (a convex combination) and takes its distance to
a third point. For static points the script times FrechetDist's `fp_` functions; for
vector points, loop versions with the same arithmetic in the same order, the
point-to-segment distance in one pass. Each kernel runs on all pairs of two random
walks of 100 vertices, in a fresh process per point type and dimension.

```sh
julia example_higher/point_types.jl            # D = 2 ... 2048, about 3 minutes
julia example_higher/point_types.jl 64,128,256 # chosen dimensions
```

Nanoseconds per call (AMD Ryzen 7 5700G, Julia 1.13.3, one thread):

| D | `dist` static | vector | `dps` static | vector | `comb` static | vector | compile static | vector |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 2 | 2.0 | 2.0 | 3.8 | 9.2 | 2.0 | 18.0 | 0.14 s | 0.20 s |
| 8 | 2.1 | 2.3 | 8.4 | 15.5 | 2.3 | 22.4 | 0.14 s | 0.20 s |
| 16 | 2.2 | 3.4 | 16.1 | 23.8 | 3.5 | 27.6 | 0.17 s | 0.20 s |
| 32 | 3.8 | 5.3 | 38.9 | 36.9 | 6.2 | 33.8 | 0.17 s | 0.19 s |
| 64 | 8.5 | 9.2 | 80.6 | 63.6 | 11.4 | 49.3 | 0.23 s | 0.20 s |
| 128 | 41.6 | 16.3 | 193.5 | 114.0 | 72.3 | 96.5 | 0.33 s | 0.20 s |
| 256 | 67.5 | 29.4 | 421.4 | 213.2 | 135.6 | 303.0 | 0.57 s | 0.20 s |
| 512 | 133.7 | 62.7 | 942.7 | 429.8 | 298.9 | 590.3 | 1.20 s | 0.20 s |
| 1024 | 299.2 | 125.2 | 2558.3 | 849.9 | 787.0 | 758.4 | 4.52 s | 0.20 s |
| 2048 | 712.8 | 236.0 | 7617.5 | 1697.2 | 1841.7 | 1624.6 | 13.49 s | 0.20 s |

The compile time is for these three kernels alone; all of `frechet_c_compute` takes
11–14 minutes at D = 6000 with static points, and about 3 s with `Vector` points at any D. Single timings in the nanosecond range vary by up to about
30% between runs; the trends do not.

- **Up to D = 16, static points win.** The point-to-segment distance is 1.5 to 2.4
  times faster, and creating a point 8 to 11 times: an `SVector` costs
  nothing to create, a `Vector` is an allocation.
- **From D = 32 to 96 it is close.** The point-to-segment distance is 1.05 to 1.5
  times faster with vectors, while creating points is still faster with
  static ones (1.6 to 5.5 times). Which wins depends on the mix of the two.
- **From D = 128 on, vectors win.** Distances are 2 to 3 times faster, and the
  point-to-segment distance 1.7 to 4.5 times. Creating points stays faster with static
  ones up to D = 512, and is even from D = 1024. Compile time also starts to matter:
  above D = 500 the static kernels alone need more than a second, and at D = 6000 the
  whole Fréchet code needs over ten minutes. On the hydraulic curves at D = 6000 (all
  vertex/edge pairs of 200 vertices), vector points ran 4.8 times faster (0.24 s
  against 1.15 s) and compiled in 0.4 s instead of 63 s.
- **The results agree.** Both point types gave identical sums for every kernel and
  dimension, except a last-bit difference for `dps` at D = 4: the squared distances
  use `@fastmath`, which lets the compiler choose the order of the additions.

The allocation in `comb` can be avoided by writing new points into reused buffers.
That would move the crossover for creating points well below D = 1024.

## HDM05

The [HDM05](https://resources.mpi-inf.mpg.de/HDM05/) database stores data in the same
ASF/AMC format, but its site returned "403 Forbidden" when these examples were written.
`mocap.jl` should read its `.amc` files as well.
