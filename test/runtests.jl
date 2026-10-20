using FrechetDist
using BasicCompGeometry
using Test

function example_01()
    P = PntSeq2F()
    Q = PntSeq2F()

    push!(P, point(0.0, 0), point(1.0, 0.0))
    push!(Q, point(0.0, 1.0), point(1.0, 1.0))
    return P, Q
end


function example_02()
    P = rand_pnt_seq(3, Float64, 100)
    Q = rand_pnt_seq(3, Float64, 80)

    return P, Q
end

#@testset "FrechetDist.jl" begin
@testset "Frechet Dist.jl" begin

    # Write your tests here.
    P, Q = example_01()

    @test FrechetDist.frechet_c_compute(P, Q).leash == 1.0

    P, Q = example_02()

    # Code generation - forcing to generate all different distance
    # calculations
    m_d = frechet_d_compute(P, Q)
    @test m_d.leash > 0.0

    m_d_r = frechet_d_r_compute(P, Q)
    @test m_d_r.leash > 0.0
    @test m_d.leash == m_d_r.leash

    m_ve = frechet_ve_r_compute(P, Q)
    @test m_ve.leash > 0.0

    m_e = frechet_c_compute(P, Q)
    @test m_e.leash > 0.0

    @test (1.000001 * m_e.leash) >= m_ve.leash

    m_a_2 = frechet_c_approx(P, Q, 2.0)

    @test (m_a_2.leash > 0.0)

    @test ((m_e.leash <= m_a_2.leash) && (m_a_2.leash <= 2.0*m_e.leash))

    d = frechet_dist_upper_bound(P, Q)
    @test m_e.leash <= d

    m_mono, f, P_r, Q_r = frechet_mono_via_refinement(P, Q, 1.01)
    @test 1.000001*m_mono.leash >= m_ve.leash

    w_P = frechet_width_approx(P)

    w_Q = frechet_width_approx(Q)
    @test w_P >= 0
    @test w_Q >= 0
end

@testset "Polygon_simplify_ext" begin
    # Always returns the simplification and the kept indices, also for an empty curve.
    Q, idx = FrechetDist.Polygon_simplify_ext(PntSeq2F(), 1.0)
    @test length(Q) == 0 && isempty(idx)
    P = PntSeq([point(0.0, 0.0), point(0.1, 0.0), point(5.0, 0.0), point(5.1, 0.0)])
    Q, idx = FrechetDist.Polygon_simplify_ext(P, 1.0)
    @test idx == [1, 3, 4] && length(Q) == 3

    Q, idx = FrechetDist.frechet_simplify_to_width(PntSeq2F(), 1.0)
    @test length(Q) == 0 && isempty(idx)
end

include("test_primitives.jl")
include("test_errors.jl")
include("test_termination.jl")
include("test_duplicates.jl")
include("test_contract.jl")
include("test_range.jl")
include("test_silence.jl")
include("test_hygiene.jl")
include("agreement.jl")

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
    for dim in (2, 3), (family, generate) in pair_families(; dim), seed = 1:seeds
        rng = MersenneTwister(hash((family, dim, seed)))
        Pv, Qv = generate(rng)
        failures = agreement_failures(Pv, Qv, rng)
        @test isempty(failures) ||
              (println("$family dim=$dim seed=$seed: ", failures); false)
    end
end
