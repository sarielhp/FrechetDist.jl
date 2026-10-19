# Computational Geometry Primitives

FrechetDist builds on [BasicCompGeometry](https://github.com/sarielhp/BasicCompGeometry.jl)
for its geometric types and operations: `Point` (an `SVector`), `Segment`,
`Line`, `PntSeq` (a polygonal curve), `BBox`, `dist`, `convex_comb`, and so on.
See that package's documentation for them.

## Floating-point primitives

The Frechet algorithms depend on a few primitives with specific floating-point
behavior: `convex_comb` snaps to an endpoint when the parameter is within
`1e-6` of it, nearest-point and segment-distance computations treat degenerate
cases explicitly, and distances use `@fastmath` accumulation. These are kept
here, with an `fp_` prefix, so that results do not change with the
BasicCompGeometry formulations.

```@autodocs
Modules = [FrechetDist]
Pages   = ["fp_geometry.jl"]
```

## Curve helpers

```@autodocs
Modules = [FrechetDist]
Pages   = ["curve_utils.jl"]
```
