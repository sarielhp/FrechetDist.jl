# Change log
    What had changed...

v3.1.0 Release Notes
====================
(November 26, 2026)

Added:
* Curves with `Vector` points of run-time length (`GPntSeq{AnyDim,Float64,Vector{Float64}}`, BasicCompGeometry 0.7): every algorithm runs on them (the polygon hierarchy, `ph_init`, takes static points only), compiling once for every dimension; above about D = 128 they are faster than static points, which take seconds to minutes to compile at D = 1000 to 6000. Matrix curves (`MatVecPntSeq`) and mixed point types are converted on entry. Results agree with static points up to rounding.

Changed:
* `Morphing{N,T,PT}`, `EventPoint{N,T,PT}` and the search contexts take the point type as a third parameter; `Morphing2F` is `Morphing{2,Float64,Point2F}`. Code using `Morphing{N,T}` still dispatches, but a field or array declared `Morphing{2,Float64}` is now abstractly typed.
* Static curves: same results bit for bit, same speed and compile time.
* Requires BasicCompGeometry 0.7 and DataStructures 0.19 (was 0.5.1/0.6 and 0.18); the priority queues use the 0.19 calls (`push!`, `popfirst!`, `first`) instead of the deprecated `enqueue!`, `dequeue!` and `peek`. Same results bit for bit.

v3.0.1 Release Notes
====================
(November 22, 2026)

Fixed (results change, see `docs/review-2026-11.md`):
* `Morphing_combine` on long curves: breakpoints on the middle curve were merged when within a relative 1e-6 of their arc-length position, 2 units at position 1e6, giving leashes up to 159 times too large. `fp_equal` now uses 1e-12 (M15). (`frechet_c_compute` also combines morphings; its results were correct in the cases tested.)
* `Polygon_simplify_ext` no longer ends with a zero-length edge when the last vertex was already kept (M14); `frechet_c_approx` keeps its leash and bounds, its morphings lose one event pair.
* `Morphing_from_prm` and `Morphing_combine` give each repeated vertex its own event (M13); before, all but the last copy had none, and their radius was 0.
* `morphing_profile`: `X` runs to the distance both points travel, `|P| + |Q|` of the morphing, matching the times it samples (before: to `|P|`) (M16).
* `Morphing_as_polygons_w_times`, `Morphing_sample_uniformly`, `morphing_profile` and `Morphing_as_function_w_times`: steps where nothing moves (or `P` does not move) take no time; they got made-up durations before, which shifted the samples.
* `frechet_c_compute` with `f_accept_approx` accepts within a factor `1 + tolerance`, not a fixed 1.00001 that ignored a smaller `tolerance`; same results with the default tolerance.
* `DTW_d_compute`: same price, bit for bit; between couplings of equal price it can return another one, and `iters` is now the number of cells.
* `frechet_width_approx`, `frechet_width_exceeds`, `frechet_offsets`, `frechet_palette` and `frechet_c_mono_approx_subcurve` reject NaN or infinite vertices with an `ArgumentError`; they returned NaN, `false` or a meaningless leash.

Faster (same results):
* `DTW_d_compute` by dynamic programming: 66x faster on random walks at n = 2000, 1/36 of the memory.
* `SweepDist_compute`: 1.9x faster on unrelated curves at n = 1000, a quarter of the allocation.
* `FEVER_compute(_range)`: 4x faster on similar curves at n = 1000 (the event arrays are no longer filled on every call).

Documentation:
* Every export has a docstring and is in the manual (new pages: Simplification, Sweep distance, Low level); the docs build fails on a missing one; the site has a `stable` link again.
* Docstrings and comments: typos and wording; `Morphing_combine` returns `v(u(.))`; `Morphing_monotonize` may return its input.

Other:
* One search loop each for VE and FEVER (`ve_search!`, `fever_search!`).
* The docs build uses Documenter 1.20 (was pinned at 1.4.1 by a stale `docs/Manifest.toml` that also held an absolute path to the package).
* The example shell scripts keep temporary files under `${TMPDIR:-/tmp}/$USER/FrechetDist/`.

v3.0.0 Release Notes
====================
(November 22, 2026)

Breaking:
* `SweepDist_compute_split` and `frechet_mono_via_refinement_ext` are renamed `SweepDist_compute_split!` and `frechet_mono_via_refinement_ext!`, and take the `out` vector first: `SweepDist_compute_split!(out, P, Q, iters, max_vertices)`, `frechet_mono_via_refinement_ext!(out, P, Q, f_snapshots, approx)`.
* Only `Float64` curves are accepted; other element types raise `ArgumentError` at the entry point (before, a `MethodError` or `InexactError` from inside).
* CSV, DataFrames, Hwloc, PrettyTables, TimerOutputs, CpuId, CPUSummary and Printf are no longer dependencies; the examples declare their own in `examples/Project.toml`.

Fixed (results change, see `docs/review-2026-11.md`):
* Sweep price of segments (`SweepDist_segs_p`, `Morphing_SweepDist_price`): closed form without absolute thresholds; it was wrong at large or small coordinate scales.
* `fp_dist_point_segment`: no cancellation at the far endpoint; only rounding noise snaps to zero, so the range functions no longer return an upper bound below the distance.
* `SweepDist_lb_compute` no longer reports a lower bound above a known morphing.
* `frechet_c_approx` and `frechet_c_compute` report a certified `lower_bound` and `ratio = leash / lower_bound`; `frechet_c_compute` no longer fails on curves far apart relative to their length.
* `DTW_d_compute`, `SweepDist_compute` and `SweepDist_lb_compute` return their objective in `sol_value`.
* Crashes on repeated event times, constant curves, 1-vertex curves and identical segments; argument validation throws `ArgumentError` instead of asserting.

Other:
* Widths, bounds and times accept any real number (for example an `Int` width).
* Faster: `fp_dist_segment_segment` no longer allocates (`SweepDist_lb_compute` 13x faster), concrete queue types, linear width scans in the simplifiers.
* `frechet_c_compute` restructured (named tuning constants, helper functions); same results.
* Uses BasicCompGeometry for points, curves, `norm`, `geom_length` and `prefix_lengths`.
* Follows the Julia standards: the lifecycle tools in `tools/`, explicit imports, one standalone module per test file, Aqua in the test suite, formatting with `always_for_in = true`.

v2.2.2 Release Notes
====================
* Repository hygiene, no behaviour change: stopped tracking scratch output, `*.jl.*.mem` allocation files and large datasets (now in `.gitignore`); added `AGENTS.md`; removed a stale GPS regression test and its fixtures.

v2.2.1 Release Notes
====================
* Internal cleanups, no behaviour change: removed dead debug code, reduced argument counts, added `test/snapshots/` (bit-for-bit comparison against a git revision) and `tools/` release scripts.

v2.2.0 Release Notes
====================
* Input contract: public entry points validate their curves and parameters and throw `ArgumentError`; library code no longer calls `exit()` or prints.
* Fixed endpoint snapping that distorted long edges and near-identical curves (leashes change by at most ~1e-16 relative).
* Fixed `SweepDist` for crossing configurations, and the penalty handling of the `*_ext` functions.
* The zero-distance clamp is relative to the input scale; coordinates above 1e150 are rejected.
* Internal: `frechet.jl` and `morphing.jl` were split into topic files (`frechet_refine.jl`, `frechet_c.jl`, `frechet_simplify.jl`, `morphing_prm.jl`).

V2.0.0 Release Notes
====================
(June 18, 2024)
* Added sweep-distance (continuous/discrete)
* Major performance improvements/bug fixes due to massive testing.
* Added `examples/prereq.jl` to install needed packages.
* Added `examples/download_data.jl` to download test data.
* Fixed github documentation generation.
* Added some documentation.


v1.0.0 Release Notes
====================
Initial release, included most of the algorithms and basic codes, was
added as an official Julia package.
