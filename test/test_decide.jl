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

@testset "the block sweep agrees with the cell sweep" begin
    FD = FrechetDist
    rng = MersenneTwister(1307)
    mismatches = 0
    for trial in 1:40
        D = (2, 3)[mod1(trial, 2)]
        step = (0.1, 1.0, 5.0)[mod1(trial, 3)]
        p = static_seq([s .* step for s in walk(rng, D, rand(rng, 10:70))])
        q = static_seq([s .* step for s in walk(rng, D, rand(rng, 10:70))])
        d = FD.frechet_c_compute_exact(p, q)
        n, m = length(p), length(q)
        for f in (1e-12, -1e-12, 1e-6, -1e-6, 0.1, -0.1), B in (2, 3, 8)
            r_sq = (d * (1 + f))^2
            s = FD.PointSpace(p.pnts, q.pnts)
            (n - 1 > B && m - 1 > B) || continue
            (FD.within(s, 1, 1, r_sq) && FD.within(s, n, m, r_sq)) || continue
            want = FD.decide_free_space(s, r_sq)
            for T in (FD.LazyBoxTree, FD.BoxTree, FD.SausageTree)
                mismatches +=
                    FD.decide_free_space_blocks(p.pnts, q.pnts, r_sq, B, T) != want
            end
        end
        mismatches += !frechet_c_decide_exact(p, q, d * (1 + 1e-9))
        mismatches += frechet_c_decide_exact(p, q, d * (1 - 1e-9))
    end
    # Vector points in ℝ^7, through the simplifications and the block sweep.
    for trial in 1:10
        a = walk(rng, 7, rand(rng, 20:60))
        b = walk(rng, 7, rand(rng, 20:60))
        d = FD.frechet_c_compute_exact(static_seq(a), static_seq(b))
        for f in (1e-9, -1e-9, 0.01, -0.01)
            want = f > 0
            for T in (FD.LazyBoxTree, FD.SausageTree)
                got = frechet_c_decide_exact(
                    vector_seq(a),
                    vector_seq(b),
                    d * (1 + f);
                    tree = T,
                )
                mismatches += got != want
            end
        end
    end
    @test mismatches == 0
end

@testset "the greedy walks and the far-vertex test are never wrong" begin
    # (the ball and sausage walks; the far-vertex test with box and sausage trees)
    FD = FrechetDist
    rng = MersenneTwister(1319)
    wrong = 0
    for trial in 1:60
        D = (2, 3)[mod1(trial, 2)]
        p = static_seq(walk(rng, D, rand(rng, 3:60))).pnts
        q = static_seq(walk(rng, D, rand(rng, 3:60))).pnts
        d = FD.frechet_c_compute_exact(PntSeq(p), PntSeq(q))
        tp, tq = FD.BoxTree(p), FD.BoxTree(q)
        sp, sq = FD.SausageTree(p), FD.SausageTree(q)
        for f in (1e-6, -1e-6, 0.01, -0.01, 0.3, -0.3)
            r = d * (1 + f)
            wrong += FD.greedy_ball_walk(p, q, tp.lens, tq.lens, r) && f < 0
            far = FD.far_vertex(p, q, tq, r^2) || FD.far_vertex(q, p, tp, r^2)
            far |= FD.far_vertex(p, q, sq, r^2) || FD.far_vertex(q, p, sp, r^2)
            wrong += FD.greedy_sausage_walk(p, q, r) && f < 0
            wrong += FD.greedy_sausage_walk(p, q, r, sp, sq) && f < 0
            wrong += far && f > 0
        end
    end
    @test wrong == 0
end

@testset "prepared curves (DecideCurve) give the same answers" begin
    FD = FrechetDist
    rng = MersenneTwister(1311)
    mismatches = 0
    for trial in 1:30
        D = (2, 3)[mod1(trial, 2)]
        a = static_seq(walk(rng, D, rand(rng, 2:40)))
        b = static_seq(walk(rng, D, trial <= 3 ? 1 : rand(rng, 2:40)))
        A, B = FD.DecideCurve(a, FD.BoxTree), FD.DecideCurve(b, FD.BoxTree)
        LA, LB = FD.DecideCurve(a), FD.DecideCurve(b)
        SA, SB = FD.DecideCurve(a, FD.SausageTree), FD.DecideCurve(b, FD.SausageTree)
        d =
            length(b) == 1 ? maximum(p -> sqrt(sum(abs2, p - b[1])), a) :
            frechet_c_compute(a, b).leash
        for f in (0.0, 1e-9, -1e-9, 0.1, -0.1, 10.0, -0.9)
            r = d * (1 + f)
            want = frechet_c_decide_exact(a, b, r)
            mismatches += frechet_c_decide_exact(A, B, r) != want
            mismatches += frechet_c_decide_exact(LA, LB, r) != want
            mismatches += frechet_c_decide_exact(SA, SB, r) != want
        end
    end
    @test mismatches == 0
    @test_throws ArgumentError frechet_c_decide_exact(
        FD.DecideCurve(static_seq([[0.0, 0.0], [1.0, 0.0]])),
        FD.DecideCurve(static_seq([[0.0, 0.0], [1.0, 0.0]])),
        -1.0,
    )
end

@testset "the default tree follows the dimension" begin
    FD = FrechetDist
    rng = MersenneTwister(1317)
    @test FD.default_tree(static_seq(walk(rng, 2, 5)).pnts) === FD.LazyBoxTree
    @test FD.default_tree(vector_seq(walk(rng, 100, 5)).pnts) === FD.SausageTree
    @test FD.DecideCurve(vector_seq(walk(rng, 100, 5))).tree isa FD.SausageTree
    @test FD.DecideCurve(static_seq(walk(rng, 3, 5))).tree isa FD.LazyBoxTree
end

@testset "LazyBoxTree has the boxes and lengths of BoxTree" begin
    FD = FrechetDist
    rng = MersenneTwister(1313)
    for n in (1, 2, 3, 5, 17, 64, 65, 300, 1000), D in (2, 3)
        p = static_seq(walk(rng, D, n)).pnts
        t, l = FD.BoxTree(p), FD.LazyBoxTree(p)
        @test FD.tree_lens(l) == t.lens
        # Read in any order: the top part and the subtrees below it are built as needed.
        for v in shuffle(rng, 1:(2*t.N-1))
            FD.node_box!(l, v)
            @test l.box[:, v] == t.box[:, v]
        end
    end
end

@testset "the sausages of SausageTree contain their pieces" begin
    FD = FrechetDist
    rng = MersenneTwister(1319)
    # the distance from x to the segment ab, in BigFloat
    function seg_dist(x, a, b)
        x, a, b = big.(x), big.(a), big.(b)
        u = b .- a
        L = sum(u .^ 2)
        s = L > 0 ? clamp(sum((x .- a) .* u) / L, 0, 1) : big(0)
        return sqrt(sum((x .- a .- s .* u) .^ 2))
    end
    line(rng, D, n) = (u = randn(rng, D); [t .* u .+ 1e-6 .* randn(rng, D) for t in 1:n])
    for n in (1, 2, 3, 9, 17, 65, 300, 1000), D in (2, 7), gen in (walk, line)
        p = vector_seq(gen(rng, D, n)).pnts
        t = FD.SausageTree(p)
        ok = true
        for v in shuffle(rng, 1:(2*t.N-1))
            FD.node_box!(t, v)
            rng_v = FD.node_range(t, v)
            isempty(rng_v) && (ok &= t.rad[1, v] < 0; continue)
            a, b = p[first(rng_v)], p[last(rng_v)+1]
            ok &=
                all(c -> seg_dist(p[c], a, b) <= t.rad[1, v], first(rng_v):(last(rng_v)+1))
        end
        @test ok
    end
end

@testset "the Fréchet radii of SausageTree bound the pieces" begin
    FD = FrechetDist
    rng = MersenneTwister(1321)
    # back and forth over a segment (rho about 0, phi about its length)
    function zigzag(rng, D, n)
        u = randn(rng, D)
        return [
            (isodd(c ÷ 9) ? 9 - c % 9 : c % 9) .* u .+ 1e-6 .* randn(rng, D) for
            c in 0:(n-1)
        ]
    end
    line(rng, D, n) = (u = randn(rng, D); [t .* u .+ 1e-6 .* randn(rng, D) for t in 1:n])
    # d_F(p[i:k], segment p[i] p[k]) <= phi, up to the tolerance of the oracle
    covers(p, i, k, phi) = oracle_continuous(p[i:k], [p[i], p[k]]) * (1 - 1e-9) <= phi
    for n in (9, 40, 70), D in (2, 7), gen in (walk, line, zigzag)
        p = [Vector(x) for x in vector_seq(gen(rng, D, n)).pnts]
        t = FD.SausageTree(p)
        ok = true
        for v in shuffle(rng, 1:(2*t.N-1))
            FD.node_box!(t, v)
            rng_v = FD.node_range(t, v)
            isempty(rng_v) && continue
            ok &= t.rad[1, v] <= t.rad[2, v]
            ok &= covers(p, first(rng_v), last(rng_v) + 1, t.rad[2, v])
        end
        for _ in 1:20  # ranges that are not nodes, as the walk reads them
            i, k = minmax(rand(rng, 1:n), rand(rng, 1:n))
            i < k || continue
            up, low = FD.frechet_bounds(t, i, k)
            ok &= covers(p, i, k, up) && low <= FD.segment_frechet(p, i, k)
            ok &= covers(p, i, k, FD.segment_frechet(p, i, k))
        end
        @test ok
    end
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
