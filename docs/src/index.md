# FrechetDist.jl
*Compute the Fréchet distance between curves.*
## Package Features
- Provides algorithm for computing the discrete, VE, regular Frechet
  distance between polygonal curves. Supports also the retractable
  version.
  

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
