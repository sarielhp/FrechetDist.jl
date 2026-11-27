# Installation

FrechetDist needs Julia 1.10 or later. In the Julia REPL:

```julia
using Pkg
Pkg.add("FrechetDist")                                           # the registered release
Pkg.add(url = "https://github.com/sarielhp/FrechetDist.jl")      # the newest version
```

The registry is updated only occasionally; the newest version, including
`frechet_c_compute_exact`, is installed from GitHub. Curves are point sequences of
[BasicCompGeometry](https://github.com/sarielhp/BasicCompGeometry.jl), which is
installed with it; add it to your environment too, to build curves:

```julia
Pkg.add("BasicCompGeometry")
```

## Threads

`frechet_c_compute_exact` builds its tables on all of Julia's threads; the other
functions run on one thread. Julia's thread count is fixed when it starts, so start it
with

```sh
julia -t auto
```

or set `JULIA_NUM_THREADS=auto` in the environment. The results do not depend on the
number of threads.

## From the source

To run the examples, clone the repository:

```sh
git clone https://github.com/sarielhp/FrechetDist.jl.git
cd FrechetDist.jl
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

The examples have their own environments (`examples/Project.toml`,
`example_higher/Project.toml`); see [Examples](examples.md) and
[High-dimensional curves: comparison with Fred](high_dimensional.md).

Next: [Quick start](quickstart.md).
