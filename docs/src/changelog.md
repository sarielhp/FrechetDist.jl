# Change log
    What had changed...

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
