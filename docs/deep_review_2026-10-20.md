# FrechetDist deep review — 2026-10-20

Scope: `src/frechet.jl`, `morphing.jl`, `fever.jl`, `SweepDist.jl`, `frechet_discrete.jl`,
`fp_geometry.jl`, `curve_utils.jl`, `palette.jl`, `polygon_hierarchy.jl`, `DistFunc.jl`,
`utils.jl`, `FrechetDist.jl` at commit `77c2cd1` (`main`).

Method: five parallel reviewers (correctness/numerics, degenerate input, architecture,
concurrency/state/limits, API/docs/tests) were seeded with 12 falsifiable hypotheses and the
independent oracle suite (`test/oracle.jl`, `test/agreement.jl`). **Every high-severity claim
below was then reproduced independently by the reviewing session** (a probe script or a
code reading), and claims that did not survive were downgraded or dropped (see section 4).
Probe scripts and prototypes are listed in section 7. Nothing in the repository was changed
by the review itself.

## 1. Executive summary

- The core algorithms are in good shape on typical data: 3,600 + 5,200 generated cases agree
  with two independent oracles; the library is thread-safe, deterministic, and does not mutate
  its inputs (0 differences in 4,480 concurrent results).
- **One silent-wrong-answer defect on realistic data (F1):** endpoint snapping at 1e-6 of an
  edge changes the matched point by up to 1e-6·L; with a long edge and short detail
  `frechet_c_compute` and `frechet_c_approx` return a wrong value *and certify it* (`ratio=1`).
  Example: edge length 1e6, bump at 0.5: result 0.583 vs exact 0.300.
- **Crashes and hangs on valid input:** `frechet_c_compute` throws an `AssertionError` on some
  inputs with duplicate consecutive vertices (F2); `frechet_simplify_to_cardin` loops forever on
  collinear curves (F3); an exported function can `exit(-1)` the whole process (F4).
- **A second silent-wrong-answer:** `Morphing_SweepDist_price` is wrong by up to 2× when the
  matched points coincide at an interior time (F5, a dead local variable).
- There is no input contract: empty, 1-vertex, NaN, Inf, `Float32`, `Int` inputs each fail in a
  different way, or silently succeed (F6). The library also prints debug text to stdout in
  several paths.
- Best payoff-to-risk: the bit-identical block (F2's cause, F3, F4, F6 validation, debug
  prints, `a == 0` branch). Highest value but behavior-changing: F1 (snapping) and F5.
- Maintainability: 12 functions over the complexity limit; a prototype shows that deleting dead
  `f_debug` code and extracting helpers takes the offenders from 12 to 3, bit-identically.

## 2. Findings (verified)

### F1 — Edge snapping makes long edges inaccurate and the certified bounds wrong
- **Severity:** high — silent wrong value, with a wrong certificate. **Confidence:** high
  (reproduced independently). **Category:** bug.
- **Location:** `src/fp_geometry.jl:47-52` (`fp_convex_comb`: snaps `t ∈ [0, 1e-6]` to `p` and
  `t ∈ (0.99999, 1]` to `q`), `src/frechet.jl:185-191` (`f_r_create_event`).
- **What's wrong:** heap values are computed analytically, but the stored event points are
  snapped. For an edge of length L the matched point moves by up to 1e-6·L (1e-5·L at the top
  end), so the morphing's leash exceeds the heap value and the "lower bound" is invalid.
- **Evidence (reproduced):** `P=[(0,0),(L,0)]`, `Q=[(0,0),(0.5,0.3),(L,0)]`; exact distance 0.300:

  | L | c_compute | c_approx (lower bound, ratio) | ve_r |
  |---|---|---|---|
  | 1e5 | 0.300 | 0.300 (ratio 1) | 0.300 |
  | 1e6 | **0.583** | **0.583 (0.583, ratio 1)** | **0.583** |
  | 1e7 | 0.583 | 0.583 | 0.583 |

  Reviewer measurements: long-edge walks with jump 1e4 fail 3/60 seeds, 1e6 7/60 (max error
  +193%); a curve against a jittered copy fails 25/80 at jitter 1e-6. The suite hides this
  because it compares with a `1e-9·scale` slack against scales around 1.
- **Why it matters:** this is the same family as the 2024 bug report (a certified value above
  the truth); anything with long edges and fine detail (sparse GPS fixes, coarse simplifications)
  can be affected.
- **Fix:** snap by distance, symmetrically (`t·|q−p| ≤ 4·eps·max|coord|`), or only at exactly
  `t == 0`/`1`; add the minimal repro as a regression test.
- **Risk:** **numerical behavior change**: needs a snapshot re-baseline and its own commit.

### F2 — `frechet_c_compute` throws on some inputs with duplicate vertices
- **Severity:** high — exception on valid input (about 1% of lattice curves with duplicates).
  **Confidence:** high (reproduced). **Category:** bug.
- **Location:** assertion `cardin(PSR) >= cardin(PS)` at `src/frechet.jl:1500`, reached from
  `:1758`.
- **Evidence:** `P=[[0,1],[-1,1],[-1,1],[-1,1]]`, `Q=[[0,1],[0,1],[-1,1],[0,3],[-1,3]]` throws in
  both `f_accept_approx` modes; 5/600 random lattice curves with duplicates threw, 0 returned a
  wrong value. `frechet_c_approx`, `frechet_ve_r_compute`, `frechet_mono_via_refinement` all
  handle the same input (oracle 2.0).
- **Root cause:** not fully traced (the refined polygons are matched back to the simplification
  by point equality, which duplicates break).
- **Fix:** carry the index map out of the refinement instead of matching points, or de-duplicate
  at entry and map indices back. **Risk:** none for inputs that work today (bit-identical).

### F3 — `frechet_simplify_to_cardin` hangs on zero-width curves and misses its target
- **Severity:** high — infinite loop on valid input. **Confidence:** high (reproduced).
- **Location:** `src/frechet.jl:1871-1911`. For a straight collinear track (or all-identical
  vertices) the width starts at 0, `R` stays at 2 vertices, and the exit test can never fire.
- **Evidence:** 60 collinear points, `sz = 10`: hang (>30 s). Contract: a 60-vertex random walk
  with `sz = 10` returns 59 vertices; a 300-vertex walk returns 190. Also an unconditional
  `println(" RRR ...")`. The function has no docstring.
- **Fix:** `m.leash == 0 && return R, R_indices`; bisect on `w` with an iteration cap; remove the
  print; document what "cardin" guarantees. **Risk:** bit-identical for curves that terminate.

### F4 — Exported functions can terminate the process
- **Severity:** high (library code must never `exit`). **Confidence:** high (reproduced: child
  exit code 255).
- **Locations:** `src/morphing.jl:721` (`check_monotone_top`, reachable from exported
  `Morphing_combine` through `parameterization_combine` with non-monotone morphings, e.g. chained
  `frechet_ve_r_compute` results), `:590` (`check_times`, reachable from the exported
  `Morphing_verify_valid` on a bad event time: verified, prints `Event time is wrong? 5.0` and
  exits). `:355` (`Morphing_recompute_leash`), `src/frechet.jl:139`, `src/fever.jl:274` are
  **not reachable from library-produced input** (reviewer analysis; the first is reachable
  only by hand-built morphings).
- **Fix:** replace each `exit(-1)` with `throw(ArgumentError(...))`, delete the dumps, document
  "monotone morphings required" on `Morphing_combine`. Bit-identical for valid input.

### F5 — `Morphing_SweepDist_price` is silently wrong (up to 2×) when matched points coincide
- **Severity:** high — silent wrong value. **Confidence:** high.
- **Location:** `src/DistFunc.jl:38-44` (`normalize` assigns `root = ...` to a dead local, so
  `f.root` stays 0 and `integral_interval` falls back to a trapezoid), `:56` (`d_eval` takes
  `sqrt` of radicands such as `-4.4e-16`, throwing `DomainError`).
- **Evidence (reproduced):** two antiparallel unit segments: `SweepDist_segs` returns 2.0, the
  numeric integral is 1.0. Reviewer: 166/200 random through-zero configurations wrong; `DomainError`
  in 3–8 of about 120 cases per producer. No test covers it.
- **Fix (reviewer-validated, 0/400 wrong afterwards):** `f.root = -f.b / (2.0 * f.a)` and
  `return v < 1e-16 ? 0.0 : sqrt(v)`. **Risk: numerical behavior change** for crossing
  configurations; separate commit.

### F6 — No input contract; inconsistent and silent failure modes
- **Severity:** medium-high. **Confidence:** high (reproduced in part).
- **Evidence:** 870 fuzz cases (15 entry points × 58 curve pairs): 669 ok, 193 exceptions
  (127 `AssertionError`, 30 `MethodError`, 25 `BoundsError`, 8 `KeyError`, 3 `ArgumentError`),
  8 hangs, 58 print to stdout. Empty or 1-vertex curves: `AssertionError`, `BoundsError`,
  `KeyError` (opaque id) or a value, depending on the function. NaN coordinates: `AssertionError`
  after printing `qr: [...]`/`t: NaN` for `c_compute`/`ve_r`, but `frechet_d_compute`/
  `DTW_d_compute` silently return a finite value that ignores the NaN (verified: 2.24 and 3.0).
  `Float32`: `MethodError` deep inside; integer `approx`: `frechet_c_approx(P, Q, 2)` throws
  `MethodError` (verified).
- **Fix:** validate at entry in every exported curve function (`ArgumentError` naming the
  function): at least 2 vertices for continuous/VE/Sweep/FEVER/mono/palette/simplify (1 for the
  discrete family), finite coordinates, magnitude ≤ 1e150, `approx > 1`, `tolerance > 0`, `n ≥ 1`;
  one `Float64` method plus a generic `ArgumentError` fallback. Bit-identical for valid input.

### F7 — Zero-length edges give NaN from the `*_range` functions and `SweepDist_lb_compute`
- **Severity:** medium. **Confidence:** high (reproduced).
- **Locations:** `fp_geometry.jl` `fp_dist_point_segment` (`a == 0` ⇒ `t = NaN` ⇒ result NaN);
  `frechet.jl:473,501`, `fever.jl:394`; `fp_dist_segment_segment`'s `@assert rk > 0`.
- **Evidence:** two curves whose last edge is zero-length: `frechet_ve_r_compute` = 0.1 (correct),
  `frechet_ve_r_compute_range` and `FEVER_compute_range` return `(NaN, NaN)`. The main
  algorithms are immune because NaN sorts last in the heap and the final leash skips NaN.
- **Fix:** an `a == 0` branch (`val = c`) in `fp_dist_point_segment` (bit-identical on 200k random
  calls with `a > 0`); `rk <= 1` in the segment-segment routine.

### F8 — Absolute constants misbehave at extreme scales
- **Severity:** medium (extreme inputs). **Confidence:** high for tiny scales, medium for huge.
- `val <= 1e-20` in `fp_dist_point_segment` clamps to zero: reproduced, coordinates of about
  1e-10 or smaller give 26–28% error (1e-300: returns 0). Scales 1e-9 … 1e153 are correct.
- Huge scales: squared distances overflow near 1e154: `c_compute` hangs at scale 2e154
  (reproduced) and throws `AssertionError` from 5e154 up. (An `Inf` coordinate throws; it does
  not hang.)
- **Fix:** short term, reject |coordinate| > 1e150 and document the minimum scale; later a
  relative clamp (`1e-20·max(a, c)`, behavior change near zero distances).

### F9 — Penalty handling in the VE machinery
- **Severity:** medium (exported `*_ext` API; no wrong result found in the default paths).
- `Morphing_recompute_leash` (`morphing.jl:347`) subtracts `pes[i].penalty` **twice**
  (`… - m.pes[i].penalty - m.pes[i].penalty`, verified by reading), never subtracts the `qes`
  penalty, and `ve_event_value` subtracts it once. The vertex-vertex branch (`frechet.jl:107`) uses
  `p_offs[i] - q_offs[j]` while both edge branches sum.
- Settled: vertex-vertex events occur only at the start and end, whose offsets are nonzero in
  `c_compute` (250 hits in 814 runs), yet flipping the sign gave bit-identical `c_compute` output
  over 7,158 runs; it matters through the exported `frechet_ve_r_compute_ext` (offsets `[0,3]`:
  leash 9.0 with `-`, 0.0 with `+`).
- **Fix:** use `+`, subtract the penalty once (and both sides consistently). Needs an input where
  the result changes as a test; classify as a behavior change.

### F10 — Contract and quality issues (medium/low)
| ID | Issue | Sev | Location |
|---|---|---|---|
| F10a | `frechet_ve_r_compute_range` throws an opaque `KeyError` and prints both curves when `upper_bound` is below the VE distance (verified: 167/168 at 0.9·VE; bound == VE exactly fails 5/168 by rounding; ≥ VE·(1+1e-9) always works). `FEVER_compute_range` falls back to a full computation instead | med | `frechet.jl:366-376, 527` |
| F10b | `tolerance` is ignored on the dominant exit: the "simplification not profitable" path calls `frechet_mono_via_refinement(P, Q)` with a hard-coded 1.00001 (604/800 calls) | med | `frechet.jl:1707` |
| F10c | `frechet_dist_upper_bound` tests `cardin(poly_a) < 2` twice (typo, never `poly_b`) and is asymmetric/invalid: for a 1-vertex first argument it returned 0.63 though the true distance is 7.07 (verified) | med | `frechet.jl:1257-1264` |
| F10d | No iteration caps in `frechet_mono_via_refinement_ext`, `frechet_c_approx`, `frechet_c_compute` (max 17 rounds seen in about 10k calls) | low | `frechet.jl:1070, 1310, 1590` |
| F10e | `frechet_c_compute` leaves `ratio = 0.0` in 300/300 results | low | `frechet.jl` |
| F10f | Debug `println` in library paths (`frechet_palette_level`, `simplify_to_cardin`, range functions, `morphing.jl:733`, `SweepDist.jl:94`) | low | several |
| F10g | `ph_approx` is a history-dependent cache (17/25 widths differ between forward and reverse call order; certificates hold) and returns the internal stored curve (aliasing) | low | `polygon_hierarchy.jl` |
| F10h | `Morphing` aliases its input curves (`m.P === P`) | low | all producers |
| F10i | `check_no_nan` checks coordinates 1 and 2 only; `EID` packs `j` into 30 bits unchecked (unreachable in memory); non-`const` global type aliases; `ph_init` requires a pointless lock argument | low | `morphing.jl:593`, `frechet.jl:19-39` |
| F10j | Quadratic-space routines allocate before validating size: measured bytes per cell `d_compute` 8, FEVER 32, `ve_r` 76–330, SweepDist about 500–600; FEVER at 20k×20k hit a bare `OutOfMemoryError` | med | `fever.jl`, `SweepDist.jl` |

### F11 — Architecture and maintainability
- Duplicated event machinery across `frechet.jl`, `fever.jl`, `SweepDist.jl`, `frechet_discrete.jl`
  (including a dead duplicate `discrete_frechet_sample`); the two VE engines break ties differently
  (equal leash in 900/900 cases, different event paths in 537/900, monotonized leash differs in
  18/750 pairs by up to 46%).
- Complexity: deleting dead `f_debug` code is bit-identical and takes `frechet_c_compute` 42 → 13,
  `frechet_c_approx` 52 → 21, `FEVER_compute_range` 19 → 15. Prototype (bit-identical on both
  snapshots, tests pass): offenders from 12 (plus depth 5) to 3. Plans for the remaining three
  (`frechet_c_approx`, `event_sequences_extract`, `add_points_along_seg`) are in the architecture
  findings file.
- `frechet.jl` (1,953 lines) and `morphing.jl` (1,123 lines) can be split into 5 and 4 files
  along existing seams (verified as pure moves, no file over 700 lines).

### F12 — API, documentation, tests
- 56 exports, 30 without a docstring; inconsistent naming (`Morphing_*`, `frechet_*`, `FEVER_*`,
  `ph_*`), return shapes (`Morphing`, tuples, `(m, f, P, Q)`) and accuracy arguments; 16 false or
  stale docstring claims (tabulated in `findings_api.md`); `docs/src/low_level.md` is not in the
  navigation.
- Coverage gaps: hierarchy/palette, `Morphing_*` utilities, mono refinement, offsets, input
  validation, thread safety. A proposed test set exists as runnable scripts (see section 7).
- The comparison snapshot `golden.jl` misses FEVER, the `*_range` functions, `frechet_d_r_compute`,
  `SweepDist_lb_compute`, `mono_via_refinement_*` and the combine paths; a supplemental snapshot
  with 1,140 entries (`golden_ext.jl`) was produced and should be adopted before the refactors.

## 3. Ordered roadmap

Order matters: stopping after any step must leave the library no worse. Steps 1–6 are
bit-identical on valid input and independently safe; steps 7–9 change numerical results and each
needs its own commit and a re-baselined snapshot.

**Now**
1. Adopt `golden_ext.jl` next to `golden.jl` as the refactor safety net (no code change).
2. `fp_dist_point_segment` `a == 0` branch; `fp_dist_segment_segment` `rk <= 1` (F7).
3. Replace the five `exit(-1)` calls with exceptions; delete the dumps (F4).
4. F3: zero-width guard and iteration cap in `frechet_simplify_to_cardin`; fix the `ph`/debug prints.
5. F2: make `frechet_c_compute` robust to duplicates (index map, not point equality).
6. F6: validation at entry (new `src/validate.jl`), plus `*_range` fallback (F10a), `poly_b` check
   (F10c), `tolerance` plumbing (F10b), iteration caps (F10d), regression tests for each.

**Next**
7. F1: snapping by distance, with the minimal repro as a test (behavior change).
8. F5: `DistFunc.jl` fix with a quadrature test (behavior change).
9. F9: penalty sign and single subtraction; F8 tiny-scale clamp (behavior changes).
10. Refactors: dead debug code deletion, then helper extractions per the prototype; split the two
    large files.

**Later**
11. API cleanup with deprecations (target about 24 exports; plan in `findings_api.md`: breaks land at
    2.1.1, 2.2.0, 3.0.0); documentation fixes; size limits for the quadratic routines.

## 4. Claims that did not survive verification, and non-issues

- **Seeded point 1 (NaN from duplicates) is largely a non-issue** for the main algorithms: NaN
  event keys sort last, cells are bypassable, and the final leash skips NaN (1,280 duplicate
  cases and 147 grid cases: 0 failures). It does leak through the `*_range` and `SweepDist_lb`
  paths (F7).
- **Thread safety, determinism, input mutation: disproved** (0/4,480 concurrent results differ;
  identical SHA-1 over three processes on 420 tie-heavy instances; 0/42 input mutations). Only
  `ph_approx`'s cache-like behavior and curve aliasing remain (F10g, F10h).
- **Downgraded:** "ve_r_range fails 575/900" is a precondition violation (bound below the VE
  distance), not a wrong result: medium; "returns 0.0 for NaN" and "frechet_dist_upper_bound
  returns 0.0" were inexact (actual behaviors above); the huge-scale hang is a narrow band
  (about 2e154) and `Inf` coordinates throw rather than hang: medium; the claimed hang at
  `1e154` did not reproduce on a different curve pair.
- **Unreachable `exit(-1)` sites:** `frechet.jl:139`, `fever.jl:274`, `morphing.jl:355` (library-
  produced input only).
- **Inputs that work** (agree with the oracle): -0.0, closed curves, identical curves, a zero-length
  edge in one curve only, 2-vertex curves, repeated sub-sequences, 3 vs 1200 vertices, interior
  edge-edge optima, ties, long straight runs, collinear back-and-forth, aspect ratio 1e6, GPS offsets
  1e7–1e9.
- Known and accepted, not reported: single-step discrete/DTW couplings; the 1e-5 tolerance in the
  `frechet_c_approx` lower bound (note that F1 shows a case where the *error is not within the
  tolerance*); `SweepDist` leash semantics.

## 5. What I would look at with more time
- A random fuzz of `Morphing_combine` with monotonized chains at larger sizes.
- Termination of `frechet_c_compute` under adversarial tolerance sequences (no cap exists).
- The palette/hierarchy extraction against a Fréchet oracle at large n.
- Whether `frechet_c_compute` and the VE engines should share one engine (F11).

## 6. Evidence note
Reproduced by the reviewing session: F1, F2, F3, F4 (both sites), F5, F7, F8, F9 (code reading),
F10a, F10c, and the NaN/short-curve behaviors of F6. Taken from the reviewers without a separate
rerun: counts such as 4,480 concurrent results, the 870-case fuzz distribution, and the
prototype measurements in F11.

## 7. Evidence artifacts
Under `/tmp/claude-10533/-home-sariel-prog-24-FrechetDist/39cac193-6e1e-45a5-8c6d-3494cf113905/scratchpad/review/`:
`prompt.md` (review prompt), `findings_api.md`, `findings_architecture.md`, my verification
probes `verify_*.jl`, the reviewers' probes (`p_snap_minimal.jl`, `p_bigjump.jl`, `p_scale.jl`,
`combine_chain.jl`, `sweep_segs.jl`, `range_nan.jl`, `fuzz2.jl`, `simplify_one.jl`, `thr.jl`,
`hier*.jl`, …), `arch/golden_ext.jl` with its baseline `base_ext.jls`, the prototypes
`arch/combined2/prototype_all.diff`, `arch/split/`, `arch/nodebug/`, and the proposed tests
`tests/t_input_contract.jl`, `t_simplify_palette_hierarchy.jl`, `t_morphing_utils.jl`, `t_threads.jl`.

## 8. Outcomes of the fix pass (branch `fix/deep-review-2026-10`)

| Finding | Status |
|---|---|
| F1 snapping | Fixed (`fp_convex_comb` snaps only at exactly 0/1; `f_r_create_event` snaps by distance relative to coordinate scale; `fp_dist_point_segment` recomputes the residual when the quadratic cancels). Behaviour change: leashes move ≤1e-16 relative; SweepDist price on D3-n120 moves ~2e-6 relative. Regression test at L up to 1e8. |
| F2, F3, F4, F6, F7 | Fixed with tests (earlier commits). |
| F5 | Fixed (`f.root` stored, negative radicands clamped). Test: antiparallel segments give 1.0 (was 2.0). Snapshots bit-identical. |
| F8 | Tiny scale: clamp made relative (`1e-20·max(a,c)`); correct to about 1e-100 (below ~1e-150 squared distances underflow). Huge scale: coordinates above 1e150 are rejected. |
| F9 | Fixed: penalty subtracted once; vertex-vertex branch uses `+`. Snapshots bit-identical (affects only exported `*_ext`). |
| F10a–d, f–j | Fixed (earlier commits). |
| F10e | Not changed: `ratio`/`lower_bound` are documented defaults of `Morphing` for exact computations. |
| F11 | Refactors: see git log. |
