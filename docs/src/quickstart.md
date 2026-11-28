# Quick start

(See [Installation](installation.md) first.)

## Two curves in the plane

```julia
using FrechetDist, BasicCompGeometry

P = PntSeq([point(0.0, 0.0), point(1.0, 1.0), point(2.0, 0.0)])
Q = PntSeq([point(0.0, 0.2), point(2.0, 0.2)])

frechet_c_compute_exact(P, Q)      # 0.8: the exact continuous Fréchet distance
frechet_c_decide_exact(P, Q, 0.8)  # true: is the distance at most 0.8?

m = frechet_c_compute(P, Q)        # a matching, to relative tolerance 1e-5
m.leash, m.lower_bound             # (0.8, 0.8): the distance lies between the two

frechet_c_approx(P, Q, 2.0).leash  # a 2-approximation, much faster on large curves
frechet_d_compute(P, Q).leash      # the discrete Fréchet distance of the vertices
DTW_d_compute(P, Q).sol_value      # dynamic time warping
```

`frechet_c_compute_exact` returns a number. Most other functions return a `Morphing`,
the matching that realizes the distance: `m.leash` is its longest leash, `m.pes` and
`m.qes` the matched points on the two curves, and, where the algorithm computes them,
`m.lower_bound` is a certified lower bound and `m.ratio = m.leash / m.lower_bound`. See
[Morphing](morphing.md) for what can be done with it (monotonizing, composing,
sampling it in time, animating it).

## High-dimensional curves

High-dimensional curves should use `Vector` points, which are compiled once for every
dimension (static points, `Point{D}`, are compiled for each D):

```julia
using FrechetDist
using BasicCompGeometry: AnyDim, GPntSeq

const Curve = GPntSeq{AnyDim,Float64,Vector{Float64}}
curve(M::Matrix{Float64}) = Curve([M[:, j] for j in axes(M, 2)])  # one column per vertex

A, B = rand(1000, 200), rand(1000, 300)   # two curves in ℝ^1000
P, Q = curve(A), curve(B)

d = frechet_c_compute_exact(P, Q)         # start Julia with -t auto to use all cores
m = frechet_c_compute(P, Q, true, 1e-9)   # tolerance 1e-9 instead of the default 1e-5
```

For timings on real curves in ℝ^963, ℝ^6000 and ℝ^19200, see [High-dimensional
curves: comparison with Fred](high_dimensional.md).

## Which function to use

| You want | Use |
|---|---|
| The exact distance, as a number | `frechet_c_compute_exact` |
| The distance to a given precision (default relative 1e-10), fast | `frechet_c_compute_bisect` |
| The distance with a matching that realizes it, or a guaranteed tolerance | `frechet_c_compute` (default tolerance 1e-5) |
| A quick estimate with guarantees | `frechet_c_approx` with `approx` = 1.1, 2, … |
| Whether the distance is at most r | `frechet_c_decide_exact` |
| Distances of point sequences (time series, samples) | `frechet_d_compute`, `DTW_d_compute` |
| A sum (integral) instead of a maximum | `SweepDist_compute_refine_mono` |
| Fewer vertices, within a given distance | `frechet_simplify_to_width`, `frechet_simplify_w_exp` |
