# Change log
    What had changed...

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
