# FrechetDist

[![Build Status](https://github.com/sarielhp/FrechetDist.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/sarielhp/FrechetDist.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/sarielhp/FrechetDist.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/sarielhp/FrechetDist.jl)
[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://sarielhp.github.io/FrechetDist.jl/dev/)

FrechetDist is a Julia library for computing the [Fréchet
distance](https://en.wikipedia.org/wiki/Fr%C3%A9chet_distance) between polygonal
curves in Euclidean space, in low dimensions (the plane, 3D) and in high dimensions
(thousands of coordinates). Informally, the Fréchet distance is the shortest leash that
lets a person and a dog walk along the two curves from start to end, both moving only
forward.

**Performance.** Approximating the distance is very fast, and `frechet_c_compute` gets
it to any chosen tolerance (10⁻⁵ by default) together with a matching that realizes it.
The exact distance, `frechet_c_compute_exact`, is not as fast: it takes quadratic time
and memory in the number of vertices. On high-dimensional curves it is still much faster
than the alternatives (below).

**Website: <https://sarielhp.github.io/FrechetDist.jl/dev/>**

## Features

- **Continuous Fréchet distance:** exact (no tolerance), to a chosen tolerance with
  the matching and a certified lower bound, approximate within any factor, and the exact
  decision procedure (is the distance at most r?).
- **High-dimensional curves, fast and exact.** On real data in ℝ^963 (freeway traffic),
  ℝ^6000 (hydraulic sensors) and ℝ^19200 (video), the exact distance of all pairs of
  each set takes 1.3 s, 160 s and 4.8 s: 85 to 143 times faster than
  [Fred](https://pypi.org/project/Fred-Frechet/), which also returns distances up to 3%
  too low for 22 of the 172 pairs.
- **Any dimension:** static points for low dimensions, `Vector` points (compiled once
  for all dimensions) for high ones.
- **Variants:** the discrete Fréchet distance, retractable versions (which minimize the
  leash at every moment), the VE-Fréchet distance, dynamic time warping, and the sweep
  distance (continuous dynamic time warping).
- **Simplification:** curves with fewer vertices within a given Fréchet distance,
  palettes and hierarchies of simplifications.
- **Morphings:** the matchings that realize the distances, which can be monotonized,
  composed, sampled in time and animated.

## Documentation

- [Installation](https://sarielhp.github.io/FrechetDist.jl/dev/installation/)
- [Quick start](https://sarielhp.github.io/FrechetDist.jl/dev/quickstart/): first
  examples, and which function to use for what
- [Implemented algorithms](https://sarielhp.github.io/FrechetDist.jl/dev/algorithms/):
  every function, grouped by problem
- [High-dimensional curves: comparison with Fred](https://sarielhp.github.io/FrechetDist.jl/dev/high_dimensional/)
- [Examples](https://sarielhp.github.io/FrechetDist.jl/dev/examples/), and
  [visual examples of output](https://sarielhp.org/p/24/frechet_ve/examples/)
- [Changelog](https://sarielhp.github.io/FrechetDist.jl/dev/changelog/)

The sources of these pages are in [`docs/src/`](docs/src/).

## Research project

This software was written as part of a research project involving Sariel Har-Peled,
Benjamin Raichel, and Eliot Robson. A paper describing the algorithms and ideas used in
this package is upcoming, and a link will be posted once it is available.
