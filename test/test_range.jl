using FrechetDist
using BasicCompGeometry
using Random
using Test

const FD = FrechetDist

include(joinpath(@__DIR__, "agreement.jl"))

@testset "Range functions with any upper bound" begin
    # frechet_ve_r_compute_range threw a KeyError (and printed both curves) when the
    # bound was below the VE distance, or equal to it up to rounding; FEVER_compute_range
    # falls back to a full computation. Both must now agree for every bound.
    for seed = 1:25
        rng = MersenneTwister(seed)
        P = mkseq(gen_walk(rng, rand(rng, 4:30)))
        Q = mkseq(gen_walk(rng, rand(rng, 4:30); step = 2.0))
        ve = FD.frechet_ve_r_compute(P, Q).leash
        for bound in (0.5 * ve, 0.9 * ve, ve, ve * (1 + 1e-9), 1.01 * ve, 2 * ve, 100 * ve)
            lo, hi = FD.frechet_ve_r_compute_range(P, Q, bound)
            @test isfinite(lo) && isfinite(hi)
            @test lo ≈ ve rtol = 1e-9
            @test hi >= lo * (1 - 1e-12)
            flo, fhi = FD.FEVER_compute_range(P, Q, bound)
            @test flo ≈ lo rtol = 1e-9
        end
    end
end

@testset "frechet_c_compute honors a tight tolerance" begin
    # The "simplification not profitable" exit used to call the monotone refinement with
    # a hard-coded 1.00001, so tolerance = 1e-12 still gave results off by up to 1e-5.
    for seed in (4, 18, 38, 39, 125, 126, 173, 194, 244, 282, 310, 347)
        rng = MersenneTwister(seed)
        Pv = gen_walk(rng, rand(rng, 4:30))
        Qv = gen_walk(rng, rand(rng, 4:30); step = 2.0)
        exact = oracle_continuous(Pv, Qv)
        leash = FD.frechet_c_compute(mkseq(Pv), mkseq(Qv), false, 1e-12).leash
        @test exact * (1 - 1e-9) <= leash <= exact * (1 + 1e-9)
    end
end
