module TestAgreement

# The Frechet algorithms against the independent oracles of oracle.jl.
using FrechetDist
using BasicCompGeometry
using Random
using Test

include(joinpath(@__DIR__, "agreement.jl"))

@testset "Oracles" begin
    # Two parallel unit segments, one at distance 1.
    P = [[0.0, 0.0], [1.0, 0.0]]
    Q = [[0.0, 1.0], [1.0, 1.0]]
    @test oracle_continuous(P, Q) ≈ 1.0
    @test oracle_discrete(P, Q) == 1.0
    # A curve against itself, and against its reverse (the endpoints are swapped).
    W = gen_walk(MersenneTwister(1), 12)
    @test oracle_continuous(W, W) == 0.0
    @test oracle_continuous(W, reverse(W)) >= oracle_dist(W[1], W[end]) - 1e-12
    # Single steps can only make the discrete distance larger.
    A, B = gen_walk(MersenneTwister(2), 15), gen_walk(MersenneTwister(3), 11)
    @test oracle_discrete(A, B; diagonal = false) >= oracle_discrete(A, B)
    @test oracle_continuous(A, B) <= oracle_discrete(A, B) * (1 + 1e-12)
end

# FrechetDist's discrete routines use single-step couplings (no diagonal moves), so they
# can exceed the textbook discrete Frechet distance and dynamic time warping price.
@testset "Discrete variants are single-step couplings" begin
    P = mkseq([[0.0, 0.0], [10.0, 0.0]])
    Q = mkseq([[0.0, 1.0], [10.0, 1.0]])
    @test FrechetDist.frechet_d_compute(P, Q).leash ≈
          oracle_discrete(vecs(P), vecs(Q); diagonal = false)
    @test oracle_discrete(vecs(P), vecs(Q); diagonal = false) >
          oracle_discrete(vecs(P), vecs(Q))
end

# The Frechet algorithms against independent oracles on random and adversarial curve
# pairs. FRECHET_STRESS_SEEDS raises the number of seeds per family (default 3).
@testset "Agreement with independent oracles" begin
    seeds = parse(Int, get(ENV, "FRECHET_STRESS_SEEDS", "3"))
    for dim in (2, 3), (family, generate) in pair_families(; dim), seed in 1:seeds
        rng = MersenneTwister(hash((family, dim, seed)))
        Pv, Qv = generate(rng)
        failures = agreement_failures(Pv, Qv, rng)
        @test isempty(failures) ||
              (println("$family dim=$dim seed=$seed: ", failures); false)
    end
end

end # module
