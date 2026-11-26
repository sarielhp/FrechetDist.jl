module TestPointTypes

# Curves with `Vector` points (`GPntSeq{AnyDim}`) and matrix curves (`MatVecPntSeq`)
# against the same curves with static points. The two paths agree up to rounding;
# where rounding changes a tie (which optimal matching is returned) the test compares
# quantities that do not depend on the tie.

using FrechetDist
using BasicCompGeometry
using Random
using StaticArrays: MVector
using Test

const FD = FrechetDist

function walk(rng, D, n)
    cur = zeros(D)
    M = zeros(D, n)
    for j in 1:n
        cur .+= rand(rng, D) .- 0.5
        M[:, j] .= cur
    end
    return M
end

static_seq(M, ::Val{D}) where {D} = PntSeq([Point{D,Float64}(c) for c in eachcol(M)])
vector_seq(M) = GPntSeq{AnyDim,Float64,Vector{Float64}}([Vector(c) for c in eachcol(M)])

const OPS = [
    "frechet_c_compute"=>(A, B, a, b) -> frechet_c_compute(A, B).leash,
    "frechet_c_approx"=>(A, B, a, b) -> frechet_c_approx(A, B, 1.01).leash,
    "frechet_ve_r_compute"=>(A, B, a, b) -> frechet_ve_r_compute(A, B).leash,
    "frechet_d_compute"=>(A, B, a, b) -> frechet_d_compute(A, B).leash,
    "frechet_d_r_compute"=>(A, B, a, b) -> frechet_d_r_compute(A, B).leash,
    "DTW_d_compute"=>(A, B, a, b) -> DTW_d_compute(A, B).sol_value,
    "SweepDist_compute"=>(A, B, a, b) -> SweepDist_compute(a, b).sol_value,
    "SweepDist_lb_compute"=>(A, B, a, b) -> SweepDist_lb_compute(a, b).leash,
    "FEVER_compute"=>(A, B, a, b) -> FEVER_compute(a, b).leash,
    "mono_via_refinement"=>(A, B, a, b) -> frechet_mono_via_refinement(a, b)[1].leash,
    "dist_upper_bound"=>(A, B, a, b) -> frechet_dist_upper_bound(A, B),
    "simplify_to_width"=>(A, B, a, b) -> length(frechet_simplify_to_width(A, 0.5)),
    "simplify_to_cardin"=>(A, B, a, b) -> length(frechet_simplify_to_cardin(A, 40)),
    "simplify_w_exp"=>(A, B, a, b) -> length(frechet_simplify_w_exp(A, 0.5)),
    "width_approx"=>(A, B, a, b) -> frechet_width_approx(A),
    "monotonize"=>(A, B, a, b) -> Morphing_monotonize(frechet_ve_r_compute(A, B)).leash,
]

@testset "Vector and matrix points agree with static points (D = $D)" for D in (2, 3, 4, 46)
    rng = MersenneTwister(5 + D)
    Ms = (walk(rng, D, 120), walk(rng, D, 110), walk(rng, D, 40), walk(rng, D, 36))
    S = map(M -> static_seq(M, Val(D)), Ms)
    V = map(vector_seq, Ms)
    C = map(MatVecPntSeq, Ms)
    for (name, f) in OPS
        s, v = f(S...), f(V...)
        @test isapprox(s, v; rtol = 1.0e-9) || (name, s, v)
        @test f(C...) == v   # a matrix curve is computed as its Vector-point copy
    end

    # The same morphing, priced with either point type.
    m = frechet_ve_r_compute(S[3], S[4])
    tovec(e) = FD.EventPoint(Vector(e.p), e.i, e.type, e.t, e.penalty)
    mv = FD.Morphing_init(V[3], V[4], map(tovec, m.pes), map(tovec, m.qes))
    @test Morphing_SweepDist_price(mv) ≈ Morphing_SweepDist_price(m) rtol = 1.0e-12
    @test Morphing_sample_uniformly(mv, 9) ≈ Morphing_sample_uniformly(m, 9) rtol = 1.0e-12
end

@testset "Vector-point kernels agree with the static ones" begin
    rng = MersenneTwister(17)
    for D in (2, 5, 46), _ in 1:200
        scale = rand(rng, (1.0e-6, 1.0, 1.0e6))
        pts = [randn(rng, D) .* scale for _ in 1:4]
        rand(rng) < 0.2 && (pts[4] = pts[3] + (pts[2] - pts[1]))  # parallel segments
        rand(rng) < 0.1 && (pts[2] = pts[1])                      # a degenerate segment
        sp = map(p -> Point{D,Float64}(p), pts)
        close(x, y) = isapprox(x, y; rtol = 1.0e-9, atol = 1.0e-12 * scale)
        @test close(
            FD.fp_dist_point_segment(pts[1:3]...),
            FD.fp_dist_point_segment(sp[1:3]...),
        )
        @test close(FD.fp_dist_segment_segment(pts...), FD.fp_dist_segment_segment(sp...))
        @test FD.fp_nn_point(pts[1:3]...) ≈ FD.fp_nn_point(sp[1:3]...)
        @test close(FD.SweepDist_segs(pts...), FD.SweepDist_segs(sp...))
    end
    p, q = [1.0, 2.0, 3.0], [4.0, 6.0, 3.0]
    @test FD.fp_convex_comb(p, q, 0.0) === p
    @test FD.fp_convex_comb(p, q, 1.0) === q
    # No allocation, except in the rank test for (nearly) parallel segments.
    r, s = [0.0, 0.0, 1.0], [1.0, 1.0, 0.0]
    # Measured in specialized functions (Julia 1.10 boxes the result at top level);
    # coverage counters allocate.
    seg_seg(a, b, c, d) = @allocated FD.fp_dist_segment_segment(a, b, c, d)
    pnt_seg(a, b, c) = @allocated FD.fp_dist_point_segment(a, b, c)
    seg_seg(p, q, r, s), pnt_seg(r, p, q)
    no_cov = Base.JLOptions().code_coverage != 0
    @test seg_seg(p, q, r, s) == 0 skip = no_cov
    @test pnt_seg(r, p, q) == 0 skip = no_cov
    @test_throws DimensionMismatch FD.fp_dist(p, [1.0, 2.0])
end

@testset "Curves the algorithms reject or convert" begin
    P3 = GPntSeq{3,Float64,Vector{Float64}}([rand(3) for _ in 1:5])
    @test_throws ArgumentError frechet_c_compute(P3, P3)
    @test_throws ArgumentError frechet_width_approx(P3)
    @test_throws ArgumentError frechet_c_compute(
        vector_seq(rand(3, 5)),
        vector_seq(rand(4, 5)),
    )
    @test_throws ArgumentError frechet_c_compute(
        static_seq(rand(2, 5), Val(2)),
        vector_seq(rand(3, 5)),
    )
    # Other point types: the kernels could not push the points they create.
    M5 = rand(3, 5)
    views = GPntSeq{AnyDim,Float64,typeof(view(M5, :, 1))}([view(M5, :, j) for j in 1:5])
    @test_throws ArgumentError frechet_c_compute(views, views)
    mvecs = GPntSeq{3,Float64,MVector{3,Float64}}([MVector{3}(rand(3)) for _ in 1:5])
    @test_throws ArgumentError frechet_c_compute(mvecs, mvecs)
    ragged = GPntSeq{AnyDim,Float64,Vector{Float64}}([rand(3) for _ in 1:4])
    push!(ragged.pnts, rand(2))   # bypasses the length check of push!
    @test_throws ArgumentError frechet_c_compute(ragged, ragged)

    # Mixed representations of the same curve: converted, not rejected.
    M = rand(3, 8)
    @test frechet_c_compute(static_seq(M, Val(3)), MatVecPntSeq(M)).leash == 0.0
    @test frechet_c_compute(static_seq(M, Val(3)), vector_seq(M)).leash == 0.0
    # The matrix is read, never written or resized.
    M0 = copy(M)
    frechet_simplify_to_width(MatVecPntSeq(M), 0.1)
    @test M == M0
    # Curves of a type the algorithms take, with wrong further arguments.
    @test_throws MethodError frechet_c_approx(vector_seq(M), vector_seq(M), "x")
    @test_throws MethodError frechet_simplify_to_width(vector_seq(M), "x")
end

end
