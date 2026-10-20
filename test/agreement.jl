# Agreement checks: FrechetDist's algorithms against the independent oracles of
# oracle.jl, on random and adversarial curve pairs. `agreement_failures` returns a
# list of human-readable failures (empty when everything agrees).

using FrechetDist
using BasicCompGeometry
using Random

const FD = FrechetDist

include(joinpath(@__DIR__, "oracle.jl"))

vecs(P) = [Vector{Float64}(p) for p in P.pnts]
mkseq(v) = PntSeq{length(v[1]),Float64}([Point{length(v[1]),Float64}(p...) for p in v])

#######################################################################
# Curve generators; all take an RNG and return a vector of coordinate vectors.

function gen_walk(rng, n; dim = 2, step = 1.0)
    cur = zeros(dim)
    pts = Vector{Vector{Float64}}()
    for _ = 1:n
        cur = cur .+ (rand(rng, dim) .- 0.5) .* step
        push!(pts, copy(cur))
    end
    return pts
end

"Dense, nearly collinear GPS-like track: long straight runs, small noise, a few turns."
function gen_track(rng, n; dim = 2, noise = 1e-3)
    pts = Vector{Vector{Float64}}()
    pos, dir = zeros(dim), normalize_dir(randn(rng, dim))
    for _ = 1:n
        rand(rng) < 0.03 && (dir = normalize_dir(dir .+ 0.8 .* randn(rng, dim)))
        pos = pos .+ dir .* (0.5 + rand(rng))
        push!(pts, pos .+ noise .* randn(rng, dim))
    end
    return pts
end

normalize_dir(v) = (n = sqrt(sum(abs2, v)); n == 0 ? [1.0; zeros(length(v) - 1)] : v ./ n)

"Points on one line, visited in random order along it (so the curve backs up)."
function gen_collinear(rng, n; dim = 2)
    dir = normalize_dir(randn(rng, dim))
    t = 0.0
    pts = Vector{Vector{Float64}}()
    for _ = 1:n
        t += rand(rng) < 0.25 ? -rand(rng) : 2rand(rng)
        push!(pts, t .* dir)
    end
    return pts
end

"Dense zigzag of fixed amplitude along a line."
function gen_zigzag(rng, n; dim = 2, amplitude = 0.3)
    return [
        vcat(
            i * 0.2,
            (isodd(i) ? amplitude : -amplitude) * (0.5 + rand(rng)),
            zeros(dim - 2),
        ) for i = 1:n
    ]
end

gen_grid(rng, n; dim = 2) = [Float64.(rand(rng, 0:3, dim)) for _ = 1:n]
gen_circle(rng, n; dim = 2) =
    [(a = 2pi * i / n; vcat(cos(a), sin(a), zeros(dim - 2))) for i = 0:(n-1)]

with_duplicates(rng, pts) = reduce(vcat, [rand(rng) < 0.3 ? [p, p] : [p] for p in pts])
jitter(rng, pts, s) = [p .+ s .* randn(rng, length(p)) for p in pts]
scaled(pts, s) = [p .* s for p in pts]

"A named family of curve pairs: `(rng) -> (P, Q)`."
function pair_families(; dim = 2)
    return [
        (
            "walks",
            rng -> (
                gen_walk(rng, rand(rng, 2:40); dim),
                gen_walk(rng, rand(rng, 2:40); dim, step = 3.0),
            ),
        ),
        (
            "near",
            rng -> (p = gen_walk(rng, rand(rng, 3:40); dim); (p, jitter(rng, p, 0.05))),
        ),
        (
            "subsample",
            rng -> (p = gen_walk(rng, rand(rng, 6:60); dim); (p, p[1:rand(rng, 2:3):end])),
        ),
        (
            "duplicates",
            rng -> (
                with_duplicates(rng, gen_walk(rng, rand(rng, 2:25); dim)),
                gen_walk(rng, rand(rng, 2:25); dim),
            ),
        ),
        ("identical", rng -> (p = gen_walk(rng, rand(rng, 2:30); dim); (p, copy(p)))),
        ("reversed", rng -> (p = gen_walk(rng, rand(rng, 2:30); dim); (p, reverse(p)))),
        ("two-vertex", rng -> (gen_walk(rng, 2; dim), gen_walk(rng, rand(rng, 2:20); dim))),
        (
            "grid",
            rng ->
                (gen_grid(rng, rand(rng, 2:25); dim), gen_grid(rng, rand(rng, 2:25); dim)),
        ),
        (
            "circles",
            rng -> (
                gen_circle(rng, rand(rng, 3:30); dim),
                jitter(rng, gen_circle(rng, rand(rng, 3:30); dim), 0.02),
            ),
        ),
        (
            "tiny",
            rng -> (
                scaled(gen_walk(rng, rand(rng, 2:25); dim), 1e-8),
                scaled(gen_walk(rng, rand(rng, 2:25); dim), 1e-8),
            ),
        ),
        (
            "huge",
            rng -> (
                scaled(gen_walk(rng, rand(rng, 2:25); dim), 1e8),
                scaled(gen_walk(rng, rand(rng, 2:25); dim), 1e8),
            ),
        ),
        (
            "tracks",
            rng -> (
                gen_track(rng, rand(rng, 20:120); dim),
                gen_track(rng, rand(rng, 20:120); dim, noise = 1e-2),
            ),
        ),
        (
            "collinear",
            rng -> (
                gen_collinear(rng, rand(rng, 3:60); dim),
                gen_collinear(rng, rand(rng, 3:60); dim),
            ),
        ),
        (
            "zigzag",
            rng -> (
                gen_zigzag(rng, rand(rng, 5:80); dim),
                gen_zigzag(rng, rand(rng, 5:80); dim),
            ),
        ),
        (
            "parallel-tracks",
            rng -> (
                p = gen_track(rng, rand(rng, 20:120); dim);
                (p, jitter(rng, [q .+ 5.0 for q in p], 1e-3))
            ),
        ),
    ]
end

#######################################################################
# Checks

# Maximum distance over the matched event pairs of a morphing: its true leash.
function morphing_leash(m)
    return maximum(oracle_dist(a.p, b.p) for (a, b) in zip(m.pes, m.qes))
end

function morphing_problems(m, name, Pv, Qv, scale)
    problems = String[]
    try
        FD.Morphing_verify_valid(m)
    catch e
        push!(problems, "$name: invalid morphing ($(first(sprint(showerror, e), 80)))")
    end
    isempty(m.pes) && return push!(problems, "$name: empty morphing")
    start_ok =
        oracle_dist(m.pes[1].p, Pv[1]) <= 1e-9 * scale &&
        oracle_dist(m.qes[1].p, Qv[1]) <= 1e-9 * scale
    end_ok =
        oracle_dist(m.pes[end].p, Pv[end]) <= 1e-9 * scale &&
        oracle_dist(m.qes[end].p, Qv[end]) <= 1e-9 * scale
    start_ok && end_ok ||
        push!(problems, "$name: morphing does not start/end at the curve endpoints")
    real = morphing_leash(m)
    abs(real - m.leash) <= 1e-9 * scale + 1e-12 * real ||
        push!(problems, "$name: reported leash $(m.leash) but its events give $real")
    return problems
end

function run_case(name, f)
    try
        return f(), nothing
    catch e
        return nothing, "$name threw $(typeof(e)): $(first(sprint(showerror, e), 110))"
    end
end

function check_discrete!(fails, P, Q, Pv, Qv, scale)
    # FrechetDist's discrete routines use single-step couplings (no diagonal moves).
    disc = oracle_discrete(Pv, Qv; diagonal = false)
    standard = oracle_discrete(Pv, Qv)
    disc >= standard * (1 - 1e-12) ||
        push!(fails, "oracle inconsistency: staircase below standard")
    m, err = run_case("frechet_d_compute", () -> FD.frechet_d_compute(P, Q))
    err === nothing || return push!(fails, err)
    isapprox(m.leash, disc; rtol = 1e-12, atol = 1e-14 * scale) ||
        push!(fails, "frechet_d_compute = $(m.leash), discrete oracle = $disc")
    append!(fails, morphing_problems(m, "frechet_d_compute", Pv, Qv, scale))
    r, err = run_case("frechet_d_r_compute", () -> FD.frechet_d_r_compute(P, Q))
    err === nothing &&
        !isapprox(r.leash, disc; rtol = 1e-12, atol = 1e-14 * scale) &&
        push!(fails, "frechet_d_r_compute = $(r.leash), discrete oracle = $disc")
    err === nothing || push!(fails, err)
end

function check_continuous!(fails, P, Q, Pv, Qv, exact, scale)
    for accept in (false, true)
        m, err =
            run_case("frechet_c_compute($accept)", () -> FD.frechet_c_compute(P, Q, accept))
        err === nothing || (push!(fails, err); continue)
        rel = exact == 0 ? m.leash / scale : (m.leash - exact) / exact
        (rel >= -1e-9 && rel <= 2e-5) || push!(
            fails,
            "frechet_c_compute($accept) = $(m.leash), exact = $exact (relative error $rel)",
        )
        append!(fails, morphing_problems(m, "frechet_c_compute($accept)", Pv, Qv, scale))
    end
    for approx in (1.0001, 1.1, 2.0)
        m, err =
            run_case("frechet_c_approx($approx)", () -> FD.frechet_c_approx(P, Q, approx))
        err === nothing || (push!(fails, err); continue)
        lower = m.leash / m.ratio
        slack = 1e-9 * scale
        # `ratio == 1` means equal up to frechet_c_compute's tolerance of 1e-5.
        (
            lower <= exact * (1 + 2e-5) + slack &&
            exact <= m.leash + slack &&
            m.leash <= approx * exact + slack
        ) || push!(
            fails,
            "frechet_c_approx($approx): lower=$lower upper=$(m.leash) ratio=$(m.ratio), exact=$exact",
        )
        append!(fails, morphing_problems(m, "frechet_c_approx($approx)", Pv, Qv, scale))
    end
end

function check_relaxations!(fails, P, Q, Pv, Qv, exact, scale)
    ve, err = run_case("frechet_ve_r_compute", () -> FD.frechet_ve_r_compute(P, Q))
    err === nothing || push!(fails, err)
    if err === nothing
        ve.leash <= exact * (1 + 1e-9) + 1e-9 * scale || push!(
            fails,
            "frechet_ve_r_compute = $(ve.leash) exceeds the exact distance $exact",
        )
        append!(fails, morphing_problems(ve, "frechet_ve_r_compute", Pv, Qv, scale))
    end
    mono, err = run_case(
        "frechet_mono_via_refinement",
        () -> FD.frechet_mono_via_refinement(P, Q, 1.01),
    )
    if err === nothing
        mono[1].leash >= exact * (1 - 1e-9) - 1e-9 * scale || push!(
            fails,
            "frechet_mono_via_refinement = $(mono[1].leash) below the exact distance $exact",
        )
    else
        push!(fails, err)
    end
    # SweepDist optimizes a sum of leash lengths and may back up along a curve, so its
    # `leash` has no fixed relation to the Frechet distance; only check its structure.
    sw, err = run_case("SweepDist_compute", () -> FD.SweepDist_compute(P, Q))
    err === nothing ?
    append!(fails, morphing_problems(sw, "SweepDist_compute", Pv, Qv, scale)) :
    push!(fails, err)
    fever, err = run_case("FEVER_compute", () -> FD.FEVER_compute(P, Q))
    if err !== nothing
        push!(fails, err)
    elseif ve !== nothing
        isapprox(fever.leash, ve.leash; rtol = 1e-9, atol = 1e-12 * scale) || push!(
            fails,
            "FEVER_compute = $(fever.leash) but frechet_ve_r_compute = $(ve.leash)",
        )
    end
end

function check_dtw!(fails, P, Q, Pv, Qv, scale)
    m, err = run_case("DTW_d_compute", () -> FD.DTW_d_compute(P, Q))
    err === nothing || return push!(fails, err)
    price = sum(oracle_dist(a.p, b.p) for (a, b) in zip(m.pes, m.qes))
    best = oracle_dtw(Pv, Qv; diagonal = false)
    isapprox(price, best; rtol = 1e-9, atol = 1e-12 * scale) || push!(
        fails,
        "DTW_d_compute price $price differs from the optimal single-step warping price $best",
    )
end

"Width approximation of random subranges, against the exact distance to the chord."
function check_width!(fails, rng, P, Pv, scale)
    length(Pv) >= 3 || return
    for _ = 1:3
        i = rand(rng, 1:(length(Pv)-2))
        j = rand(rng, (i+2):length(Pv))
        w, err = run_case("frechet_width_approx", () -> FD.frechet_width_approx(P, i:j))
        err === nothing || (push!(fails, err); continue)
        chord = oracle_continuous(Pv[i:j], [Pv[i], Pv[j]])
        (
            w >= chord * (1 - 1e-9) - 1e-9 * scale &&
            w <= 2 * chord * (1 + 1e-9) + 1e-9 * scale
        ) || push!(
            fails,
            "frechet_width_approx($i:$j) = $w, exact distance to the chord = $chord",
        )
    end
end

"Simplifications must stay within their requested width of the original curve."
function check_simplifications!(fails, rng, P, Pv, scale)
    length(Pv) >= 3 || return
    total = sum(oracle_dist(Pv[i], Pv[i+1]) for i = 1:(length(Pv)-1))
    for frac in (0.02, 0.1, 0.4)
        w = frac * total / max(1, length(Pv) - 1) * 3
        w > 0 || continue
        r, err =
            run_case("frechet_simplify_to_width", () -> FD.frechet_simplify_to_width(P, w))
        err === nothing || (push!(fails, err); continue)
        isempty(r[1].pnts) && continue
        d = oracle_continuous(Pv, vecs(r[1]))
        d <= w * (1 + 1e-9) + 1e-9 * scale || push!(
            fails,
            "frechet_simplify_to_width($w): result is at distance $d from the curve",
        )
        plt = FD.frechet_palette(P)
        q, idx, wq = FD.frechet_palette_approx(P, plt, w)
        dq = oracle_continuous(Pv, vecs(q))
        dq <= max(w, wq) * (1 + 1e-9) + 1e-9 * scale || push!(
            fails,
            "frechet_palette_approx($w): result is at distance $dq (reported width $wq)",
        )
    end
end

"""
    agreement_failures(Pv, Qv, rng)

Run every checkable property of the Frechet algorithms on the curve pair and return
the list of failures.
"""
function agreement_failures(Pv, Qv, rng)
    fails = String[]
    P, Q = mkseq(Pv), mkseq(Qv)
    exact = oracle_continuous(Pv, Qv)
    scale = max(exact, maximum(oracle_dist(p, Pv[1]) for p in Pv), 1e-300)
    check_discrete!(fails, P, Q, Pv, Qv, scale)
    check_continuous!(fails, P, Q, Pv, Qv, exact, scale)
    check_relaxations!(fails, P, Q, Pv, Qv, exact, scale)
    check_dtw!(fails, P, Q, Pv, Qv, scale)
    check_width!(fails, rng, P, Pv, scale)
    check_simplifications!(fails, rng, P, Pv, scale)
    return fails
end
