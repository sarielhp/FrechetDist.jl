# FrechetDist.jl
*Compute the Fréchet distance between curves.*
## Package Features
- Provides algorithm for computing the discrete, VE, regular Frechet
  distance between polygonal curves. Supports also the retractable
  version.

## Point types and high dimensions

Curves are `BasicCompGeometry` point sequences. Every algorithm takes either kind:

| Curve | Points | Use |
|---|---|---|
| `PntSeq{D,Float64}` | static (`Point{D,Float64}`) | low dimensions (fastest up to D = 16) |
| `GPntSeq{AnyDim,Float64,Vector{Float64}}` | `Vector{Float64}` | high dimensions (faster from D = 128; compiles once for every D) |
| `MatVecPntSeq(M)` | the columns of a matrix | converted on entry, see below |

Both curves of a call must have the same dimension and coordinate type (`Float64`).
Results with `Vector` points agree with static points up to rounding, not bit for bit:
a last-bit difference can select another optimal morphing, with the same leash.

A `MatVecPntSeq` (or a static curve paired with a `Vector`-point curve) is copied once,
at the entry, into `Vector` points, since the algorithms create new points; the matrix
is never written. Morphings returned for it have `Vector` points. Results share `Vector` points with
the input curves (static points are copied by value): writing into a point,
`P[i][k] = x`, changes it in both. BasicCompGeometry's `setindex!`, `translate!` and
`scale!` replace points and are safe. The polygon hierarchy (`ph_init`) takes static
points only. A sequence with
`Vector` points must have the dimension parameter `AnyDim`:
`GPntSeq{3,Float64,Vector{Float64}}` is rejected with an `ArgumentError` that says how
to rewrap it without copying. Other point types (views, `MVector`s) are rejected too:
the algorithms create points and push them into the curves.

`Morphing{N,T}` and `EventPoint{N,T}` take the point type as a third parameter
(`Morphing{2,Float64,Point2F}`, `Morphing2F`). The two-parameter forms still dispatch,
but are no longer concrete types: a field or array declared with them is abstractly
typed.

## Function Documentation
```@docs 
frechet_c_compute
frechet_c_approx
frechet_dist_upper_bound
FrechetDist.frechet_width_approx
```

## VE-Frechet distance and refinement

The VE-Frechet distance considers only vertex-edge events. Its retractable version is
a lower bound on the Frechet distance, and monotonizing its morphing gives an upper
bound; refining the curves along the morphing closes the gap.

```@docs
frechet_ve_r_compute
frechet_ve_r_compute_ext
frechet_ve_r_compute_range
FEVER_compute
FEVER_compute_range
frechet_mono_via_refinement
frechet_mono_via_refinement_ext!
frechet_mono_via_refinement_delta
frechet_refinement
```
