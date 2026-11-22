# Simplification

Simplifying a curve replaces it by one with fewer vertices, within a given Frechet
distance of it.

## Simplification to a width

```@docs
frechet_simplify_to_width
frechet_simplify_w_exp
frechet_simplify_to_cardin
simplify_morphing_sensitive
```

## Palettes and hierarchies

A palette is a sequence of simplifications of a curve at decreasing widths; a polygon
hierarchy stores them so that an approximation of any accuracy can be read off.

```@docs
frechet_palette
frechet_palette_level
frechet_palette_approx
PolygonHierarchy
ph_init
ph_approx
FrechetDist.frechet_measure
```
