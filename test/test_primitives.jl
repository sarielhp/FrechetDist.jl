using FrechetDist
using BasicCompGeometry
using Test

const FD = FrechetDist

@testset "Zero-length segments in the distance primitives" begin
    p = point(1.0, 2.0)
    q = point(4.0, 6.0)
    # The distance to a zero-length segment is the distance to its single point.
    @test FD.fp_dist_point_segment(q, p, p) ≈ 5.0
    @test FD.fp_dist_point_segment(p, p, p) == 0.0
    @test FD.fp_dist_segment_segment(p, p, q, q) ≈ 5.0
    @test FD.fp_dist_segment_segment(p, q, q, q) == 0.0
    @test FD.fp_dist_segment_segment(q, q, p, q) == 0.0
    @test FD.fp_dist_segment_segment(p, p, p, p) == 0.0
end

@testset "Range functions with zero-length edges" begin
    # Both last edges have zero length; the VE distance is 0.1.
    P = PntSeq([point(0.0, 0.0), point(1.0, 0.0), point(1.0, 0.0)])
    Q = PntSeq([point(0.0, 0.1), point(1.0, 0.1), point(1.0, 0.1)])
    ve = FD.frechet_ve_r_compute(P, Q).leash
    @test ve ≈ 0.1
    for range in
        (FD.frechet_ve_r_compute_range(P, Q, 5.0), FD.FEVER_compute_range(P, Q, 5.0))
        lo, hi = range
        @test isfinite(lo) && isfinite(hi)
        @test lo ≈ ve
        @test hi >= lo
    end
end
