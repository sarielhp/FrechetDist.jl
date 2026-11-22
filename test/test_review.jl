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
