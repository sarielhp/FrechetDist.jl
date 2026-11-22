# Morphing

A *morphing* is an encoding of a "matching" between two curves. It is
the output of the various Frechet distance computation functions.  The
package uses the geometric types of BasicCompGeometry (`Point`, `Segment`, `PntSeq`).

## Types
```@docs
Morphing
EventPoint
FPointType
```

## Functions
```@docs
Morphing_monotonize
Morphing_empty
Morphing_verify_valid
Morphing_combine
Morphing_extract_prm
Morphing_as_polygons
Morphing_extract_vertex_radii
Morphing_extract_offsets
Morphing_monotone_leash
Morphing_get_max_edges_err
```

## Morphings as functions of time

```@docs
FrechetDist.Morphing_as_polygons_w_times
Morphing_as_function_w_times
polygons_get_loc_at_time
Morphing_sample_uniformly
morphing_profile
```
