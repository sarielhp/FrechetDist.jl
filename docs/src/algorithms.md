# Implemented algorithms

Every exported function, grouped by problem. Each has a docstring with the details
(`?frechet_c_compute` in the REPL), also in the reference pages of this manual.
Here n and m are the numbers of vertices of the two curves, N the larger of them, and
d the dimension.

## Continuous Fréchet distance

| Function | What it computes |
|---|---|
| `frechet_c_compute_exact(P, Q)` | **The exact distance**, with no tolerance: the randomized algorithm of Har-Peled and Raichel, [*The Fréchet distance revisited and extended*](https://arxiv.org/abs/1202.5610) (Section 6). The distance is the smallest critical value of Alt and Godau (vertex–vertex, vertex–segment and monotonicity events) that the decision procedure accepts. A random sample of 4N² critical values, binary searched with the decider, leaves an interval with O(N log N) critical values with high probability; a kinetic sweep over a balanced search tree lists them, and a second binary search finds the distance. O(nmd) time to tabulate vertex–segment distances (on all threads), then O(nm) per decision. The result does not depend on the random choices or the number of threads. Returns a `Float64`. |
| `frechet_c_compute_bisect(P, Q; rtol, atol)` | **The distance to a given precision**: a bisection with the decision procedure, as FrechetLight::calcDistance of Bringmann, Künnemann and Nusser, until the interval holding the distance is at most `max(atol, rtol * d)` long (default `rtol = 1e-10`). The first decision is at the distance of the end points plus the tolerance, which settles the many pairs whose distance is that of their end points. About log2 of the ratio of the starting interval to the tolerance decisions, most of them settled by the filters; on the real curves of the benchmark of Conradi, van der Hoog and Rotenberg (SoCG 2026), faster than the C++ codes of Bringmann et al. and of CGAL at their precision. Returns a `Float64`. |
| `frechet_c_decide_exact(P, Q, r)` | **The decision procedure**: is the distance at most `r`? The free-space reachability test of Alt and Godau, O(nmd) time and O(n+m) memory, with exact shortcuts (bounds, greedy walks, a far-vertex test, simplifications, a recursive decomposition of the free space into boxes, with the boundary rules of Bringmann, Künnemann and Nusser: a box whose top or right edge is free in one interval is settled from its boundary) that settle most queries much faster. |
| `frechet_c_compute(P, Q, f_accept_approx = true, tol = 1e-5)` | **The distance to a relative tolerance, with a matching.** Simplifies both curves at a scale set by a rough approximation, computes the VE-retractable distance of the simplifications (a certified lower bound, using vertex offsets) and a monotone matching (an upper bound), and refines until the two agree to `tol`. Returns a monotone `Morphing`: `m.leash` is the leash of an actual matching, so never below the distance, and `m.lower_bound` is certified. Usually the fastest way to the distance on large curves. |
| `frechet_c_approx(P, Q, approx)` | **An approximation** within a factor `approx` (which may be larger than 2, for a rough one). Returns a monotone morphing with `m.leash ≤ approx · distance` and the certified lower bound `m.lower_bound`. |
| `frechet_mono_via_refinement(P, Q, approx)` | The monotone distance by refinement alone: computes the VE-retractable morphing, monotonizes it, and adds vertices along the morphing until the monotone leash is within `approx` of the VE-retractable distance. No simplification, so slow on huge curves. |
| `frechet_mono_via_refinement_ext!(out, P, Q, f_snapshots, approx)` | The same, recording the morphing of every refinement round (for animations). |
| `frechet_mono_via_refinement_delta(P, Q, delta, f_fix_prm)` | The same with an additive gap `delta`; optionally maps the morphing back to the original curves. |
| `frechet_refinement(P, Q, m)` | Refines both curves along a morphing: adds the points where its events lie on the edges. |
| `frechet_dist_upper_bound(P, Q)` | A rough, fast upper bound, used as a starting point by the algorithms above. |

## VE-retractable Fréchet distance

The VE-Fréchet distance considers only vertex–edge events; its *retractable* version
lets the walk move backwards and minimizes the leash at every moment, not only the
largest one. It is a lower bound on the Fréchet distance, and monotonizing its morphing
gives an upper bound. These are the building blocks of `frechet_c_compute` and
`frechet_c_approx`.

| Function | What it computes |
|---|---|
| `frechet_ve_r_compute(P, Q)` | The VE-retractable distance and its (possibly non-monotone) morphing, by a Dijkstra-like search over vertex–edge events. |
| `frechet_ve_r_compute_ext(P, Q, p_offs, q_offs, f_use_offsets)` | The same with per-vertex offsets subtracted from the event distances: used to certify lower bounds for simplified curves. |
| `frechet_ve_r_compute_range(P, Q)` | A range `(lo, hi)`: the VE-retractable distance and the leash of its monotonized morphing, a lower and an upper bound on the Fréchet distance. |
| `FEVER_compute(P, Q)`, `FEVER_compute_range(P, Q, upper_bound)` | The same distance and range with quadratic-space lookup tables instead of hashing: faster for small curves (up to about 200 vertices). |
| `frechet_c_mono_approx_subcurve(P, Q, p_indices)` | The distance between a curve and a subcurve given by vertex indices (as produced by simplification). |

## Discrete Fréchet distance and dynamic time warping

| Function | What it computes |
|---|---|
| `frechet_d_compute(P, Q)` | The discrete Fréchet distance of the two vertex sequences, by dynamic programming, O(nm). Couplings advance one index at a time (no diagonal steps), so the result can exceed the textbook discrete distance. |
| `frechet_d_r_compute(P, Q)` | The retractable discrete Fréchet distance: a Prim/Dijkstra search that minimizes the bottleneck leash at every moment, not only the largest one. |
| `frechet_d_r_compute_sample(P, Q, n, f_lopt = true)` | The continuous distance approximated discretely: samples about `n` points on the two curves (in proportion to their lengths, plus all vertices) and computes the retractable (or, with `f_lopt = false`, the standard) discrete distance. |
| `DTW_d_compute(P, Q)` | Dynamic time warping: the minimum, over discrete couplings, of the sum of the leash lengths, by dynamic programming, O(nm). The price is `m.sol_value`. |

## Sweep distance (continuous dynamic time warping)

The sweep distance of two curves is the smallest integral of the leash length over a
morphing between them.

| Function | What it computes |
|---|---|
| `SweepDist_compute(P, Q)` | A morphing minimizing the approximate sweep price, by a Dijkstra search over vertex–edge events (not necessarily monotone). |
| `SweepDist_compute_refine_mono(P, Q)` | A monotone morphing for the sweep distance: refines the curves along the morphing and recomputes until monotonizing it costs less than 1% of the curves' length. |
| `SweepDist_compute_split!(out, P, Q, iters)` | Morphings at increasing resolution, splitting every edge at its midpoint in each round. |
| `SweepDist_lb_compute(P, Q)` | A lower bound on the sweep distance. |
| `Morphing_SweepDist_price(m)`, `Morphing_SweepDist_approx_price(m)` | The sweep price of a given morphing, exactly and as a fast estimate. |

## Simplification

| Function | What it computes |
|---|---|
| `frechet_simplify_to_width(P, w)` | Greedy simplification with Fréchet error about `w`: repeatedly jumps to the farthest vertex whose shortcut stays within width `w` (binary search). Returns the curve and the indices of the kept vertices. |
| `frechet_simplify_w_exp(P, w)` | The same with an exponential search before each binary search, so a jump costs time proportional to its length. Faster when many vertices are kept. |
| `frechet_simplify_to_cardin(P, sz)` | A simplification with at least about `sz` vertices, by decreasing the allowed error. |
| `simplify_morphing_sensitive(m, factor)` | Simplifies both curves of a morphing where it has slack, keeping the vertices that determine the leash. |
| `frechet_width_approx(P, rng)` | A 2-approximation of the Fréchet distance between a subcurve `P[i:j]` and the segment `P[i] P[j]`, by a greedy monotone matching. |
| `frechet_offsets(P, p_indices)` | For every edge of a subcurve, a fast upper estimate of its distance to the part of `P` it replaces. |
| `frechet_palette(P)`, `frechet_palette_level(P, plt, level)`, `frechet_palette_approx(P, plt, w)` | A palette: a per-vertex encoding of a whole sequence of simplifications at decreasing widths, from which a simplification of any width, or of a given level, is read off in time proportional to its size. |
| `PolygonHierarchy`, `ph_init`, `ph_compute_hierarchy`, `ph_approx`, `ph_push!`, `ph_print` | A thread-safe hierarchy of cached simplifications of a curve, from which an approximation of any accuracy is extracted. |

## Morphings

Most functions return a `Morphing`: the matched event points on the two curves (`pes`,
`qes`, of type `EventPoint`), the leash, and, where the algorithm computes them, a
certified lower bound, the approximation ratio and a price.

| Function | What it does |
|---|---|
| `Morphing_monotonize(m)`, `Morphing_monotone_leash(m)` | Turns a morphing into a monotone one (never going back), or computes its leash without building it. |
| `Morphing_combine(u, v)` | The composition of two morphings, from `u.P` to `v.Q`. |
| `Morphing_verify_valid(m)`, `Morphing_empty` | Structural validation; an empty morphing. |
| `Morphing_as_polygons(m)`, `Morphing_extract_prm(m)` | The morphing as two polygons whose edges are matched directly, or as a parameterization: a monotone piecewise-linear curve in the plane. |
| `Morphing_as_function_w_times(m)`, `polygons_get_loc_at_time`, `Morphing_sample_uniformly(m, n)`, `morphing_profile(m, n)` | The morphing as a function of time: the two locations at any time, uniform samples, and the leash length along the walk. |
| `Morphing_extract_vertex_radii(m)`, `Morphing_extract_offsets(m)`, `Morphing_get_max_edges_err(m)` | Per-vertex leash bounds and the loosest pair of matched edges, used to guide refinement and simplification. |
