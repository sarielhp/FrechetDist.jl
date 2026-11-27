module TestDecide

# frechet_c_decide_exact against the independent free-space oracle of oracle.jl, against
# frechet_c_compute, on exact boundary cases, and across point types.

using FrechetDist
using BasicCompGeometry
using Random
using Test

include(joinpath(@__DIR__, "oracle.jl"))

function walk(rng, D, n)
    cur = zeros(D)
    return [copy(cur .+= rand(rng, D) .- 0.5) for _ in 1:n]
end

static_seq(v) = PntSeq([Point{length(v[1]),Float64}(p) for p in v])
vector_seq(v) = GPntSeq{AnyDim,Float64,Vector{Float64}}([Vector(p) for p in v])

const FACTORS = (0.5, 0.999, 1 - 1e-9, 1 + 1e-9, 1.001, 2.0)

@testset "frechet_c_decide_exact agrees with the oracle around the distance" begin
    rng = MersenneTwister(1301)
    disagreements = 0
    for trial in 1:60
        D = (2, 3, 7)[mod1(trial, 3)]
        a = walk(rng, D, rand(rng, 2:25))
        b = walk(rng, D, rand(rng, 2:25))
        d = oracle_continuous(a, b)
        for f in FACTORS
            r = d * f
            want = oracle_decide(a, b, r)
            disagreements += frechet_c_decide_exact(static_seq(a), static_seq(b), r) != want
            disagreements += frechet_c_decide_exact(vector_seq(a), vector_seq(b), r) != want
        end
    end
    @test disagreements == 0
end

@testset "frechet_c_decide_exact brackets frechet_c_compute" begin
    rng = MersenneTwister(1302)
    for trial in 1:30
        D = (2, 3, 40)[mod1(trial, 3)]
        a, b = walk(rng, D, rand(rng, 3:40)), walk(rng, D, rand(rng, 3:40))
        P, Q = D == 40 ? (vector_seq(a), vector_seq(b)) : (static_seq(a), static_seq(b))
        m = frechet_c_compute(P, Q, true, 1e-10)
        @test frechet_c_decide_exact(P, Q, m.leash * (1 + 1e-8))
        @test !frechet_c_decide_exact(P, Q, m.leash * (1 - 1e-6))
    end
end

@testset "frechet_c_decide_exact at exact distances" begin
    seq(pts...) = PntSeq([Point2F(p...) for p in pts])
    # Parallel segments at distance 1.
    P, Q = seq((0, 0), (1, 0)), seq((0, 1), (1, 1))
    @test frechet_c_decide_exact(P, Q, 1.0)
    @test !frechet_c_decide_exact(P, Q, prevfloat(1.0))
    # A vertex at distance 1 from the middle of a segment.
    P, Q = seq((0, 0), (2, 0)), seq((0, 0), (1, 1), (2, 0))
    @test frechet_c_decide_exact(P, Q, 1.0)
    @test !frechet_c_decide_exact(P, Q, prevfloat(1.0))
    # Backtracking: Q goes to 3 and back to 1, so P must wait; the distance is 1. Here
    # two free intervals are compared, (3 - r)/4 against (1 + r)/4, and at 1 ulp below
    # 1 both round to 0.5 (the oracle answers the same); 1e-12 below is decided.
    P, Q = seq((0, 0), (4, 0)), seq((0, 0), (3, 0), (1, 0), (4, 0))
    @test frechet_c_decide_exact(P, Q, 1.0)
    @test !frechet_c_decide_exact(P, Q, 1.0 - 1e-12)
    # Identical curves, repeated vertices, one-vertex curves.
    P = seq((0, 0), (1, 2), (1, 2), (3, 1))
    @test frechet_c_decide_exact(P, P, 0.0)
    @test frechet_c_decide_exact(seq((0, 0), (1, 0)), seq((0.5, 0)), 0.5)
    @test !frechet_c_decide_exact(seq((0, 0), (1, 0)), seq((0.5, 0)), prevfloat(0.5))
    @test frechet_c_decide_exact(seq((1, 1)), seq((1, 1)), 0.0)
    # The endpoints decide alone.
    @test !frechet_c_decide_exact(seq((0, 0), (1, 0)), seq((0, 0), (1, 3)), 2.9)
end

@testset "frechet_c_decide_exact input handling" begin
    P = PntSeq([Point2F(0, 0), Point2F(1, 0)])
    Q = PntSeq([Point2F(0, 1), Point2F(1, 1)])
    @test_throws ArgumentError frechet_c_decide_exact(P, Q, -1.0)
    @test_throws ArgumentError frechet_c_decide_exact(P, Q, NaN)
    @test_throws ArgumentError frechet_c_decide_exact(P, Q, Inf)
    @test_throws ArgumentError frechet_c_decide_exact(P, PntSeq{2,Float64}(), 1.0)
    R = PntSeq([Point{3,Float64}(0, 1, 0), Point{3,Float64}(1, 1, 0)])
    @test_throws ArgumentError frechet_c_decide_exact(P, R, 1.0)
    # Matrix curves are converted on entry.
    M = MatVecPntSeq([0.0 1.0; 0.0 0.0])
    N = MatVecPntSeq([0.0 1.0; 1.0 1.0])
    @test frechet_c_decide_exact(M, N, 1.0)
    @test !frechet_c_decide_exact(M, N, 0.99)
end

end
