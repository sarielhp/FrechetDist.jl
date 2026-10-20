using FrechetDist
using BasicCompGeometry
using Random
using Test

const FD = FrechetDist

include(joinpath(@__DIR__, "oracle.jl"))

curve(coords...) = PntSeq([point(c...) for c in coords])
vecs(P) = [Vector{Float64}(p) for p in P.pnts]

@testset "frechet_c_compute with duplicate consecutive vertices" begin
    # Minimal case found by the deep review: used to throw
    # AssertionError(cardin(PSR) >= cardin(PS)) in get_refinement_map.
    P = curve((0.0, 1.0), (-1.0, 1.0), (-1.0, 1.0), (-1.0, 1.0))
    Q = curve((0.0, 1.0), (0.0, 1.0), (-1.0, 1.0), (0.0, 3.0), (-1.0, 3.0))
    exact = oracle_continuous(vecs(P), vecs(Q))
    @test exact ≈ 2.0
    for accept in (false, true)
        m = FD.frechet_c_compute(P, Q, accept)
        @test exact * (1 - 1e-9) <= m.leash <= exact * (1 + 2e-5)
    end

    # Random lattice curves with repeated vertices (6 of 600 used to throw).
    function lattice(rng, n)
        pts = [Float64.(rand(rng, 0:3, 2)) for _ = 1:n]
        return reduce(vcat, [rand(rng) < 0.3 ? [p, p] : [p] for p in pts])
    end
    for seed = 1:300
        rng = MersenneTwister(seed)
        A, B = lattice(rng, rand(rng, 3:14)), lattice(rng, rand(rng, 3:14))
        exact = oracle_continuous(A, B)
        leash = FD.frechet_c_compute(
            PntSeq([point(p...) for p in A]),
            PntSeq([point(p...) for p in B]),
        ).leash
        @test exact * (1 - 1e-9) - 1e-12 <= leash <= exact * (1 + 2e-5) + 1e-12
    end

    # A curve of identical vertices still works against a real curve.
    same = curve((1.0, 1.0), (1.0, 1.0), (1.0, 1.0))
    other = curve((0.0, 0.0), (2.0, 0.0))
    @test FD.frechet_c_compute(same, other).leash ≈
          oracle_continuous(vecs(same), vecs(other)) rtol = 2e-5
end
