# FrechetDist.jl
*Compute the Fréchet distance between curves.*
## Package Features
- Provides algorithm for computing the discrete, VE, regular Frechet 
  distance between curves. Supports also the retractable version.
  
## Basic geometric types

```@autodocs
Modules = [FrechetDist.cg]
Order   = [:type, :function]
```

## Function Documentation
```@docs 
frechet_c_compute
frechet_c_approx
frechet_dist_upper_bound
```
