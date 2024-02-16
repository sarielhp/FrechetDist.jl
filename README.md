# FrechetDist

[![Build Status](https://github.com/YourUserNameOnGithub/FrechetDist.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/YourUserNameOnGithub/FrechetDist.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/YourUserNameOnGithub/FrechetDist.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/YourUserNameOnGithub/FrechetDist.jl)

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

[Documentation](build/docs/)

## Examples

To generate examples, you can run:

> julia examples/generate_examples.jl

This generates movies/etc in the subdirectory output/01/, output/02/,...

## Research project

This software was written as part of a research project involving
Sariel Har-Peled, Benjamin Raichel, and Eliot Robson. A paper
describing the algorithms/ideas used in this package is upcoming, and
a link would be posted once it is available.
