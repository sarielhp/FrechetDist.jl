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
    for seed = 1:60
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
