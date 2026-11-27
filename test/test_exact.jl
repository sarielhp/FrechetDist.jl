module TestExact

# frechet_c_compute_exact: the kinetic sweep against listing every critical value, the
# distance against the independent oracle of oracle.jl and frechet_c_compute, exact
# values, and the Las Vegas property (the result does not depend on the random sample).

using FrechetDist
using BasicCompGeometry
using Random
using Test

const FD = FrechetDist

include(joinpath(@__DIR__, "oracle.jl"))

walk(rng, D, n) = (cur = zeros(D); [copy(cur .+= rand(rng, D) .- 0.5) for _ in 1:n])
# Small integer coordinates: many equal distances, collinear vertices, coincidences.
grid(rng, D, n) = [Float64.(rand(rng, -2:2, D)) for _ in 1:n]
# Collinear curves that go back and forth.
line(rng, n) = [[Float64(rand(rng, 0:6)), 0.0] for _ in 1:n]

static_seq(v) = PntSeq([Point{length(v[1]),Float64}(p) for p in v])
vector_seq(v) = GPntSeq{AnyDim,Float64,Vector{Float64}}([Vector(p) for p in v])

function curve_pair(rng, trial)
    n, m = rand(rng, 2:14), rand(rng, 2:14)
    kind = trial % 4
    kind == 0 && return line(rng, n), line(rng, m)
    kind == 3 && return walk(rng, 4, n), walk(rng, 4, m)
    return grid(rng, 1 + kind, n), grid(rng, 1 + kind, m)
end

vertices(v) = FD.remove_consecutive_duplicates(static_seq(v)).pnts

@testset "kinetic sweep lists the same critical values as brute force" begin
    rng = MersenneTwister(1401)
    differ = 0
    for trial in 1:400
        a, b = curve_pair(rng, trial)
        p, q = vertices(a), vertices(b)
        lo, hi = sort(rand(rng, 2) .* 4)
        # The sweep starts from the order at lo, so values within an ulp or two of lo
        # may be missing; the search lists from below its lower end for that reason.
        above(S) = filter(>(nextfloat(lo, 8)), S)
        S1 = above(FD.exact_critical_values(p, q, lo, hi))
        S2 = above(FD.exact_critical_values(p, q, lo, hi; brute = true))
        differ += S1 != S2
    end
    @test differ == 0
end

@testset "frechet_c_compute_exact agrees with the oracle; same result for every rng" begin
    rng = MersenneTwister(1402)
    worst = 0.0
    seed_dependent = 0
    point_type_dependent = 0
    for trial in 1:200
        a, b = curve_pair(rng, trial)
        d = frechet_c_compute_exact(
            static_seq(a),
            static_seq(b);
            rng = MersenneTwister(trial),
        )
        o = oracle_continuous(a, b; rtol = 1e-15)
        worst = max(worst, abs(d - o) / max(o, floatmin()))
        d2 = frechet_c_compute_exact(
            static_seq(a),
            static_seq(b);
            rng = MersenneTwister(10^6 + trial),
        )
        seed_dependent += d2 !== d
        point_type_dependent += frechet_c_compute_exact(vector_seq(a), vector_seq(b)) !== d
    end
    @test worst <= 1e-12
    @test seed_dependent == 0
    @test point_type_dependent == 0
end

@testset "frechet_c_compute_exact against frechet_c_compute" begin
    rng = MersenneTwister(1403)
    for trial in 1:12
        D = (2, 3, 40)[mod1(trial, 3)]
        a, b = walk(rng, D, rand(rng, 3:30)), walk(rng, D, rand(rng, 3:30))
        P, Q = D == 40 ? (vector_seq(a), vector_seq(b)) : (static_seq(a), static_seq(b))
        m = frechet_c_compute(P, Q, true, 1e-10)
        d = frechet_c_compute_exact(P, Q)
        @test m.lower_bound * (1 - 1e-9) <= d <= m.leash * (1 + 1e-9)
    end
end

@testset "frechet_c_compute_exact at exact distances" begin
    seq(pts...) = PntSeq([Point2F(p...) for p in pts])
    # A vertex-vertex, a vertex-segment and a monotonicity value, all exactly 1.
    @test frechet_c_compute_exact(seq((0, 0), (1, 0)), seq((0, 1), (1, 1))) == 1.0
    @test frechet_c_compute_exact(seq((0, 0), (2, 0)), seq((0, 0), (1, 1), (2, 0))) == 1.0
    @test frechet_c_compute_exact(
        seq((0, 0), (4, 0)),
        seq((0, 0), (3, 0), (1, 0), (4, 0)),
    ) == 1.0
    # Identical curves, repeated vertices, one-vertex curves.
    P = seq((0, 0), (1, 2), (1, 2), (3, 1))
    @test frechet_c_compute_exact(P, P) == 0.0
    @test frechet_c_compute_exact(seq((0, 0), (1, 0), (0, 0)), seq((0.5, 0))) == 0.5
    @test frechet_c_compute_exact(seq((3, 4)), seq((0, 0))) == 5.0
end

@testset "frechet_c_compute_exact input handling" begin
    P = PntSeq([Point2F(0, 0), Point2F(1, 0)])
    R = PntSeq([Point{3,Float64}(0, 1, 0), Point{3,Float64}(1, 1, 0)])
    @test_throws ArgumentError frechet_c_compute_exact(P, R)
    @test_throws ArgumentError frechet_c_compute_exact(P, PntSeq{2,Float64}())
    M = MatVecPntSeq([0.0 1.0; 0.0 0.0])
    N = MatVecPntSeq([0.0 1.0; 1.0 1.0])
    @test frechet_c_compute_exact(M, N) == 1.0
end

@testset "the brute-force fallback was never needed" begin
    @test FD.EXACT_FALLBACKS[] == 0
end

end
