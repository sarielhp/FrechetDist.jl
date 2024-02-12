#  Discrete Fréchet distance
## Description
   The discrete Frechet distance dells with the easier variant, where
   the input is two sequence of points, and the distance is the
   minimum max-leash that is needed, if two "frogs" jumps in a
   synchronized fashion along the two sequences, starting at the start
   of the sequences, and end in the end of the sequences. At each
   point only one frog jumps, and it can jump only forward by one
   position in the sequence of points it is on.

  The basic solution to this is dynamic programming, which is similar
   in nature to edit-distance.

  We also present the retractable variant, which tries to use a leash
   that is as short as possible at any point in time.

  Finally, we also present a variant for the "continuous" variant,
   where one samples points along two polygons, and then compute their
   discrete Frechet distance. This is the "standard" way of using the
   discrete Frechet distance for the continuous problem. It provides
   inferior results than the continuous variants, but it is simple
   enough to implement, and it is thus provided.

## Functions documentation
```@docs
frechet_d_compute
frechet_d_r_compute
frechet_d_r_compute_sample
```
