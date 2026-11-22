# FrechetDist.jl / Low level details
## Miscellaneous

- We cover here a few low level internal functions. This is probably
only of interest if you are trying to understand the source coe of the
package.


## Functions documentation
```@docs
FrechetDist.extract_param_inner
FrechetDist.extract_refined_polygon
FrechetDist.frechet_c_mono_approx_subcurve
FrechetDist.frechet_ve_r_mono_compute
FrechetDist.Morphing_from_prm
FrechetDist.add_points_along_seg!
FrechetDist.EID
FrechetDist.SweepDist_segs
FrechetDist.frechet_offsets
FrechetDist.SweepDist_segs_p
FrechetDist.SweepDist_lb_compute
FrechetDist.frechet_width_exceeds
FrechetDist.max_leash
FrechetDist.get_point
FrechetDist.sqrt_quad_integral
FrechetDist.PrevMap
FrechetDist.PrevLookup
FrechetDist.PREVMAP_MIN_DENSE
FrechetDist.check_events
FrechetDist.event_close
```

## Argument validation

```@autodocs
Modules = [FrechetDist]
Pages   = ["validate.jl"]
```
