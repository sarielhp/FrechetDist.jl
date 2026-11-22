module TestReview

# Regression tests for the defects found by the code review of 2026-11.
using FrechetDist
using BasicCompGeometry
using Random
using Test

const FD = FrechetDist

# The value of `f()`, or the exception it throws: lets a crash fail a test
# (rather than error it) when compared against the expected value.
outcome(f) =
    try
        f()
    catch e
        e
    end

mkc(coords) = PntSeq([point(c...) for c in coords])

@testset "polygons_get_loc_at_time with repeated times" begin
    P = mkc([(0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (3.0, 0.0)])
    Q = mkc([(0.0, 1.0), (1.0, 1.0), (2.0, 1.0), (3.0, 1.0)])
    times = [0.0, 0.5, 0.5, 1.0]
    @test outcome(() -> FD.polygons_get_loc_at_time(P, Q, times, 0.25)) ==
          (point(0.5, 0.0), point(0.5, 1.0))
    @test outcome(() -> FD.polygons_get_loc_at_time(P, Q, times, 0.5)) ==
          (point(1.0, 0.0), point(1.0, 1.0))
    @test outcome(() -> FD.polygons_get_loc_at_time(P, Q, times, 0.75)) ==
          (point(2.5, 0.0), point(2.5, 1.0))

    # Continuous morphings of these walks have repeated event times.
    fails = 0
    for seed in 1:60
        rng = MersenneTwister(seed)
        A = PntSeq([point(x...) for x in eachcol(cumsum(randn(rng, 3, 30), dims = 2))])
        B = PntSeq([point(x...) for x in eachcol(cumsum(randn(rng, 3, 25), dims = 2))])
        r = outcome(() -> FD.morphing_profile(FD.frechet_c_compute(A, B), 200))
        r isa Exception && (fails += 1)
    end
    @test fails == 0

    # Degenerate inputs: a constant curve pair, and a single sample.
    C = mkc([(1.0, 1.0), (1.0, 1.0)])
    X, Y = FD.morphing_profile(FD.frechet_ve_r_compute(C, C), 5)
    @test all(iszero, Y)
    m = FD.frechet_ve_r_compute(P, Q)
    @test outcome(() -> FD.Morphing_sample_uniformly(m, 1)) isa ArgumentError
    @test size(FD.Morphing_sample_uniformly(m, 3)) == (3, 4)
end

@testset "DTW and sweep distance report their objective in sol_value" begin
    rng = MersenneTwister(5)
    A = PntSeq([point(x...) for x in eachcol(cumsum(randn(rng, 2, 30), dims = 2))])
    B = PntSeq([point(x...) for x in eachcol(cumsum(randn(rng, 2, 25), dims = 2))])
    m = FD.DTW_d_compute(A, B)
    price = sum(FD.fp_dist(a.p, b.p) for (a, b) in zip(m.pes, m.qes))
    @test m.sol_value ≈ price rtol = 1e-12
    s = FD.SweepDist_compute(A, B)
    @test s.sol_value ≈ FD.Morphing_SweepDist_approx_price(s) rtol = 1e-9
end

@testset "Degenerate but valid inputs do not crash" begin
    A = mkc([(1.0, 1.0), (1.0, 1.0)])          # a curve of length zero
    B = mkc([(0.0, 0.0), (2.0, 0.0), (2.0, 2.0)])
    for lopt in (true, false)
        @test outcome(() -> FD.frechet_d_r_compute_sample(A, B, 10, lopt).leash) ≈ sqrt(2.0)
        @test outcome(() -> FD.frechet_d_r_compute_sample(B, A, 10, lopt).leash) ≈ sqrt(2.0)
    end

    # A tiny delta used to overflow round(Int64, length / delta).
    P = mkc([(0.0, 0.0), (2.0, 0.0)])
    Q = mkc([(0.0, 0.0), (1.0, 0.0), (0.5, 0.0), (2.0, 0.0)])
    for delta in (1e-18, 1e-20)
        @test outcome(() -> FD.frechet_mono_via_refinement_delta(P, Q, delta)[1].leash) ≈
              0.25
    end

    # One-vertex curves are accepted by the simplifiers.
    one = mkc([(3.0, 4.0)])
    for f in (FD.frechet_simplify_to_width, FD.frechet_simplify_w_exp)
        r = outcome(() -> f(one, 1.0))
        @test r isa Tuple && length(r[1]) == 1 && r[2] == [1]
    end

    # Identical segments have sweep price zero (used to be NaN, then an error).
    for p in (point(0.0, 0.0), point(1.0, 2.0, 3.0))
        q = p + p .+ 1.0
        @test outcome(() -> FD.SweepDist_segs(p, q, p, q)) == 0.0
    end
    m = FD.frechet_c_compute(B, B)
    @test outcome(() -> FD.Morphing_SweepDist_price(m)) == 0.0
end

@testset "Invalid arguments raise ArgumentError" begin
    P = mkc([(0.0, 0.0), (1.0, 1.0), (2.0, 0.0)])
    Q = mkc([(0.0, 1.0), (1.0, 2.0), (2.0, 1.0)])
    AE = ArgumentError
    for bad in (NaN, -1.0, Inf)
        offs = [0.0, bad, 0.0]
        @test outcome(() -> FD.frechet_ve_r_compute_ext(P, Q, offs, zeros(3), true)) isa AE
    end
    R = mkc([(1.0, 1.0), (2.0, 0.0)])
    @test outcome(() -> FD.frechet_c_mono_approx_subcurve(P, R, [2, 3])) isa AE
    @test outcome(() -> FD.frechet_c_mono_approx_subcurve(P, R, [1, 2])) isa AE
    @test outcome(() -> FD.frechet_width_approx(P, 2:7)) isa AE
    plt = FD.frechet_palette(P)
    @test outcome(() -> FD.frechet_palette_approx(P, plt[1:2], 1.0)) isa AE
    @test outcome(() -> FD.frechet_palette_level(P, plt[1:2], 1)) isa AE
    m = FD.frechet_ve_r_compute(P, Q)
    @test outcome(() -> FD.Morphing_combine(FD.Morphing_empty(P, Q), m)) isa AE
    @test outcome(() -> FD.simplify_morphing_sensitive(m, 0.0)) isa AE
    @test outcome(() -> FD.frechet_simplify_to_cardin(P, -5)) isa AE

    # Wrong extra arguments are a MethodError, not a claimed dimension mismatch.
    @test outcome(() -> FD.frechet_c_compute(P, Q, 1e-3)) isa MethodError
    @test outcome(() -> FD.frechet_c_approx(P, Q)) isa MethodError
    P3 = mkc([(0.0, 0.0, 0.0), (1.0, 1.0, 1.0)])
    @test outcome(() -> FD.frechet_c_compute(P, P3)) isa AE

    # Event ids must still decode when j == n_q + 1 for the largest allowed n_q.
    j = FD.MAX_VERTICES + 1
    @test FD.EID_j(FD.EID(3, false, j, true)) == j
    @test FD.EID_i(FD.EID(3, false, j, true)) == 3
end

@testset "Morphing_verify_valid checks the structure of a morphing" begin
    rng = MersenneTwister(9)
    A = PntSeq([point(x...) for x in eachcol(cumsum(randn(rng, 2, 30), dims = 2))])
    B = PntSeq([point(x...) for x in eachcol(cumsum(randn(rng, 2, 25), dims = 2))])
    m = FD.frechet_ve_r_compute(A, B)
    @test outcome(() -> FD.Morphing_verify_valid(m)) === nothing
    AE = ArgumentError

    bad = deepcopy(m)
    pop!(bad.pes)                                          # length mismatch
    @test outcome(() -> FD.Morphing_verify_valid(bad)) isa AE

    bad = deepcopy(m)
    e = bad.pes[2]
    bad.pes[2] = FD.EventPoint(e.p, 999, e.type, e.t, e.penalty)  # index out of range
    @test outcome(() -> FD.Morphing_verify_valid(bad)) isa AE

    bad = deepcopy(m)
    k = findfirst(e -> e.type == FD.PT_ON_EDGE, bad.pes)
    e = bad.pes[k]
    bad.pes[k] = FD.EventPoint(e.p .+ 1.0, e.i, e.type, e.t, e.penalty)  # moved point
    @test outcome(() -> FD.Morphing_verify_valid(bad)) isa AE

    bad = deepcopy(m)
    reverse!(bad.pes)                                      # does not start at P[1]
    reverse!(bad.qes)
    @test outcome(() -> FD.Morphing_verify_valid(bad)) isa AE
end

@testset "fp_dist_segment_segment: no allocation, same values as the LAPACK solve" begin
    a, b, c, d = point(0.0, 0.0), point(1.0, 0.2), point(0.3, 1.0), point(0.5, -1.0)
    FD.fp_dist_segment_segment(a, b, c, d)
    @test (@allocated FD.fp_dist_segment_segment(a, b, c, d)) == 0

    ref =
        isdefined(FD, :fp_dist_segment_segment_lapack) ? FD.fp_dist_segment_segment_lapack :
        FD.fp_dist_segment_segment
    rng = MersenneTwister(21)
    mism = 0
    for k in 1:20_000, D in (2, 3)
        sc = 2.0^rand(rng, -60:60)
        pts = [Point{D,Float64}(randn(rng, D) .* sc) for _ in 1:4]
        if k % 5 == 0                       # lattice: exact ties, parallel, shared points
            pts = [Point{D,Float64}(rand(rng, -2:2, D)) for _ in 1:4]
        end
        r1 = FD.fp_dist_segment_segment(pts...)
        r2 = ref(pts...)
        r1 === r2 || (mism += 1)
    end
    @test mism == 0
end

@testset "Search contexts have concrete field types" begin
    for T in (FD.SDContext{2,Float64}, FD.FRContext{2,Float64})
        @test all(isconcretetype, fieldtypes(T))
    end
end

@testset "frechet_width_exceeds agrees with frechet_width_approx" begin
    rng = MersenneTwister(13)
    A = PntSeq([point(x...) for x in eachcol(cumsum(randn(rng, 2, 200), dims = 2))])
    bad = 0
    for _ in 1:3000
        i = rand(rng, 1:200)
        j = rand(rng, i:200)
        w = FD.frechet_width_approx(A, i:j)
        for v in (w, prevfloat(w), nextfloat(w), 0.5w, 2w, 0.0)
            v >= 0 || continue
            FD.frechet_width_exceeds(A, i:j, v) == (w > v) || (bad += 1)
        end
    end
    @test bad == 0
end

@testset "VE search: dense and Dict predecessor maps give identical results" begin
    conv(m) = (
        m.leash,
        m.iters,
        [(e.p, e.i, e.type, e.t) for e in m.pes],
        [(e.p, e.i, e.type, e.t) for e in m.qes],
    )
    rng = MersenneTwister(17)
    same = true
    for k in 1:40
        D = 2 + k % 2
        A = PntSeq([
            Point{D,Float64}(x) for
            x in eachcol(cumsum(randn(rng, D, rand(rng, 3:40)), dims = 2))
        ])
        B = PntSeq([
            Point{D,Float64}(x) for
            x in eachcol(cumsum(randn(rng, D, rand(rng, 3:40)), dims = 2))
        ])
        a = conv(FD.frechet_ve_r_compute(A, B))
        r = FD.frechet_ve_r_compute_range(A, B, 0.01)
        old = FD.PREVMAP_MIN_DENSE[]
        FD.PREVMAP_MIN_DENSE[] = 0
        try
            same &= (a == conv(FD.frechet_ve_r_compute(A, B)))
            same &= (r == FD.frechet_ve_r_compute_range(A, B, 0.01))
        finally
            FD.PREVMAP_MIN_DENSE[] = old
        end
    end
    @test same
end

@testset "SweepDist_lb_compute is a lower bound" begin
    # A zero dp value used to be overwritten by a larger one (dp > 0 as "set").
    P = mkc([(1.0, 3.0), (2.0, 0.0)])
    Q = mkc([(2.0, 2.0), (2.0, 0.0), (3.0, 0.0), (3.0, 0.0)])
    price = FD.Morphing_SweepDist_price(FD.SweepDist_compute_refine_mono(P, Q))
    @test FD.SweepDist_lb_compute(P, Q).sol_value <= price * (1 + 1e-12)

    rng = MersenneTwister(4)
    viol = 0
    for _ in 1:300
        A = mkc([(rand(rng, 0:3) * 1.0, rand(rng, 0:3) * 1.0) for _ in 1:rand(rng, 2:6)])
        B = mkc([(rand(rng, 0:3) * 1.0, rand(rng, 0:3) * 1.0) for _ in 1:rand(rng, 2:6)])
        lb = FD.SweepDist_lb_compute(A, B).sol_value
        pr = FD.Morphing_SweepDist_price(FD.SweepDist_compute_refine_mono(A, B))
        lb <= pr * (1 + 1e-9) + 1e-12 || (viol += 1)
    end
    @test viol == 0
end

@testset "fp_dist_point_segment beyond the far end of a long segment" begin
    s_p, s_q = point(0.0, 0.0), point(1.0e6, 0.0)
    qr = point(1.0e6 + 2.0^-10, 2.0^-10)   # exactly representable
    @test FD.fp_dist_point_segment(qr, s_p, s_q) ≈ sqrt(2.0) * 2.0^-10 rtol = 1e-12
end

@testset "Distances far below the segment length are not snapped to zero" begin
    for L in (1.0e8, 1.0e10, 1.0e12)
        P = mkc([(0.0, 0.0), (L / 2, 0.5), (L, 0.0)])
        Q = mkc([(0.0, 0.0), (L, 0.0)])
        @test FD.fp_dist_point_segment(P[2], Q[1], Q[2]) ≈ 0.5
        lo, hi = FD.frechet_ve_r_compute_range(P, Q, 10.0)
        @test hi >= 0.5 * (1 - 1e-9)
        lo, hi = FD.FEVER_compute_range(P, Q, 10.0)
        @test hi >= 0.5 * (1 - 1e-9)
    end
    # Points on a segment still give exactly zero.
    rng = MersenneTwister(2)
    nonzero = 0
    for _ in 1:10_000
        sc = 10.0^rand(rng, -6:8)
        a, b =
            point(randn(rng) * sc, randn(rng) * sc), point(randn(rng) * sc, randn(rng) * sc)
        t = rand(rng, 0:16) / 16                        # exact on-segment points
        p = a + t * (b - a)
        FD.fp_dist_point_segment(p, a, b) == 0.0 || (nonzero += 1)
    end
    @test nonzero == 0
end

@testset "Sweep price: accurate at every scale, no spurious errors" begin
    # Near-translation: the leash stays about 1 long (was classified as |linear|, got 1.0).
    p_a, p_b = point(0.0, 0.0), point(1.0, 0.0)
    q_a, q_b = point(0.0, 1.0), point(1.0 + 1.2e-7, 1.0 - 7.2e-15)
    @test FD.SweepDist_segs(p_a, p_b, q_a, q_b) ≈ 2.0 rtol = 1e-6

    # Parallel 1 m segments 100 m apart at UTM-like coordinates (was 0: fp_equal).
    for O in (0.0, 1.0e6, 4.0e6)
        @test FD.SweepDist_segs(
            point(O, O),
            point(O + 1, O),
            point(O, O + 100),
            point(O + 1, O + 100),
        ) ≈ 200.0 rtol = 1e-9
    end

    # Scale invariance: price(s x) = s^2 price(x).
    pts = (point(0.3, -0.2), point(1.1, 0.9), point(-0.4, 0.5), point(0.8, -1.3))
    base = FD.SweepDist_segs(pts...)
    for s in (1.0e-9, 1.0e-6, 1.0e4, 1.0e8)
        @test FD.SweepDist_segs((s .* p for p in pts)...) / s^2 ≈ base rtol = 1e-12
    end

    # Random walks at large coordinates used to throw DomainError in sqrt.
    fails = 0
    for seed in 1:40
        rng = MersenneTwister(seed)
        A = PntSeq([
            point(x...) for x in eachcol(1.0e4 .* cumsum(randn(rng, 2, 20), dims = 2))
        ])
        B = PntSeq([
            point(x...) for x in eachcol(1.0e4 .* cumsum(randn(rng, 2, 20), dims = 2))
        ])
        r = outcome(() -> FD.Morphing_SweepDist_price(FD.frechet_ve_r_compute(A, B)))
        (r isa Real && isfinite(r)) || (fails += 1)
    end
    @test fails == 0
end

include(joinpath(@__DIR__, "agreement.jl"))

@testset "frechet_c_approx / frechet_c_compute report certified bounds" begin
    # The shortcut path of frechet_c_approx reported lower_bound = leash, above the
    # exact distance (e.g. subsample, seed 15: 0.333226589 > 0.333225697), the main
    # path a value other than leash / ratio, and frechet_c_compute left ratio = 0.
    bad = String[]
    for (fname, gen) in pair_families(; dim = 2), seed in 1:15
        Pv, Qv = gen(MersenneTwister(seed))
        P, Q = mkseq(Pv), mkseq(Qv)
        exact = oracle_continuous(Pv, Qv)
        slack = 1e-9 * max(exact, 1e-300)
        ms = [("approx($ap)", FD.frechet_c_approx(P, Q, ap)) for ap in (1.0001, 1.1, 2.0)]
        append!(
            ms,
            [("compute($acc)", FD.frechet_c_compute(P, Q, acc)) for acc in (true, false)],
        )
        for (nm, m) in ms
            ok =
                m.lower_bound <= exact + slack &&
                exact <= m.leash + slack &&
                (exact == 0 || m.lower_bound > 0) &&
                (
                    m.lower_bound == 0 ||
                    isapprox(m.ratio, m.leash / m.lower_bound; rtol = 1e-12)
                )
            ok || push!(
                bad,
                "$fname/$seed $nm: lb=$(m.lower_bound) leash=$(m.leash) ratio=$(m.ratio) exact=$exact",
            )
        end
    end
    @test isempty(bad)
    isempty(bad) || println(join(first(bad, 5), "\n"))
end

@testset "frechet_c_compute on curves far apart relative to their length" begin
    # approx = 1 + length / (100 leash) rounded to exactly 1.0 and was rejected.
    for sh in (1.0e15, 1.0e16)
        rng = MersenneTwister(1)
        A = PntSeq([point(x...) for x in eachcol(cumsum(randn(rng, 2, 5), dims = 2))])
        B = PntSeq([
            point((x .+ [sh, 0.0])...) for x in eachcol(cumsum(randn(rng, 2, 5), dims = 2))
        ])
        r = outcome(() -> FD.frechet_c_compute(A, B).leash)
        @test r isa Real && abs(r - sh) <= 1.0e-3 * sh
    end
end

end # module
