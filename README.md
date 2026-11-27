# FrechetDist

[![Build Status](https://github.com/YourUserNameOnGithub/FrechetDist.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/YourUserNameOnGithub/FrechetDist.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/YourUserNameOnGithub/FrechetDist.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/YourUserNameOnGithub/FrechetDist.jl)

**High-dimensional curves, fast and exact.** FrechetDist handles curves in any
dimension (with `Vector` points, compiled once for all dimensions). On real data in
ℝ^963 (freeway traffic), ℝ^6000 (hydraulic sensors) and ℝ^19200 (video), its exact
continuous Fréchet distance computes all pairs of each set in 1.3 s, 160 s and 4.8 s:
85 to 143 times faster than [Fred](https://pypi.org/project/Fred-Frechet/), which
also returns distances up to 3% too low for 22 of the 172 pairs. See [High-dimensional
curves: comparison with Fred](https://sarielhp.github.io/FrechetDist.jl/dev/high_dimensional/)
([source](docs/src/high_dimensional.md)).

This package implements in Julia algorithms for computing the Frechet distance
between curves. The [Fréchet
distance](https://en.wikipedia.org/wiki/Fr%C3%A9chet_distance) is
informally the person-dog distance. Specifically, given two curves,
imagine the dog walking along one of the curves, and the person
walking on the other. Throughout the walk, there is a leash keeping
them connected, and the problem is to plan the motion that minimizes
the longest leash necessary to keep them together. For an animation
demonstrating how such a synchronized walk might look like, see
[here](https://www.youtube.com/watch?v=j5pPVC-mw_I).

This library is written in Julia and implements a number of different
algorithms for the Fréchet distance, including the "easy" discrete
variant, a retractable variant, a continuous monotone variant, etc.

[Documentation](https://sarielhp.github.io/FrechetDist.jl/dev/)

[Visual examples of output](  https://sarielhp.org/p/24/frechet_ve/examples/)

## Examples

To generate examples, you can run:

> julia examples/generate_examples.jl

This generates movies/etc in the subdirectory output/01/, output/02/,...

## Research project

This software was written as part of a research project involving
Sariel Har-Peled, Benjamin Raichel, and Eliot Robson. A paper
describing the algorithms/ideas used in this package is upcoming, and
a link would be posted once it is available.
