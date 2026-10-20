using FrechetDist
using BasicCompGeometry
using Test

const FD = FrechetDist

mk(coords) = PntSeq([point(c...) for c in coords])
two = mk([(0.0, 0.0), (1.0, 0.0)])
three = mk([(0.0, 1.0), (1.0, 1.5), (2.0, 1.0)])
one = mk([(0.0, 0.0)])
empty = PntSeq2F()

# name => (function of two curves, minimum number of vertices)
entry_points = Dict(
    "frechet_d_compute" => (FD.frechet_d_compute, 1),
    "frechet_d_r_compute" => (FD.frechet_d_r_compute, 1),
    "DTW_d_compute" => (FD.DTW_d_compute, 1),
    "frechet_ve_r_compute" => (FD.frechet_ve_r_compute, 2),
    "FEVER_compute" => (FD.FEVER_compute, 2),
    "frechet_c_compute" => (FD.frechet_c_compute, 2),
    "frechet_c_approx" => ((P, Q) -> FD.frechet_c_approx(P, Q, 1.5), 2),
    "SweepDist_compute" => (FD.SweepDist_compute, 2),
    "SweepDist_lb_compute" => (FD.SweepDist_lb_compute, 2),
    "frechet_mono_via_refinement" => (FD.frechet_mono_via_refinement, 2),
    "frechet_dist_upper_bound" => (FD.frechet_dist_upper_bound, 2),
)

function outcome(f, P, Q)
    try
        f(P, Q)
        return :value
    catch e
        return e isa ArgumentError ? :ArgumentError : Symbol(typeof(e))
    end
end

@testset "Input contract of the public entry points" begin
    for (name, (f, min_vertices)) in sort(collect(entry_points); by = first)
        @testset "$name" begin
            @test outcome(f, empty, three) == :ArgumentError
            @test outcome(f, three, empty) == :ArgumentError
            @test outcome(f, one, three) == (min_vertices == 1 ? :value : :ArgumentError)
            @test outcome(f, mk([(0.0, 0.0), (NaN, 1.0)]), three) == :ArgumentError
            @test outcome(f, three, mk([(0.0, 0.0), (Inf, 1.0)])) == :ArgumentError
            @test outcome(f, mk([(0.0, 0.0), (1.0e151, 1.0)]), three) == :ArgumentError
            # Valid input still works.
            @test outcome(f, two, three) == :value
        end
    end

    @testset "Element type and dimension" begin
        float32 = PntSeq([Point{2,Float32}(0.0f0, 0.0f0), Point{2,Float32}(1.0f0, 1.0f0)])
        three_d = PntSeq([point(0.0, 0.0, 0.0), point(1.0, 1.0, 1.0)])
        for (name, (f, _)) in sort(collect(entry_points); by = first)
            @test outcome(f, float32, float32) == :ArgumentError
            @test outcome(f, float32, three) == :ArgumentError
            @test outcome(f, three_d, three) == :ArgumentError
        end
    end

    @testset "Parameters" begin
        # An integer accuracy is fine; the accuracy must exceed 1 and be finite.
        @test FD.frechet_c_approx(two, three, 2).leash > 0
        @test_throws ArgumentError FD.frechet_c_approx(two, three, 1.0)
        @test_throws ArgumentError FD.frechet_c_approx(two, three, 0.5)
        @test_throws ArgumentError FD.frechet_c_approx(two, three, NaN)
        @test_throws ArgumentError FD.frechet_mono_via_refinement(two, three, NaN)
        @test FD.frechet_mono_via_refinement(two, three, 2)[1].leash > 0
        @test_throws ArgumentError FD.frechet_c_compute(two, three, true, -1.0)
        @test_throws ArgumentError FD.frechet_c_compute(two, three, true, NaN)
        @test FD.frechet_c_compute(two, three, true, 1e-8).leash > 0
        @test_throws ArgumentError FD.frechet_ve_r_compute_range(two, three, -1.0)
        @test_throws ArgumentError FD.FEVER_compute_range(two, three, NaN)
        @test_throws ArgumentError FD.frechet_simplify_to_width(three, -1.0)
        @test_throws ArgumentError FD.frechet_simplify_w_exp(three, NaN)
        @test_throws ArgumentError FD.frechet_ve_r_compute_ext(
            two,
            three,
            [0.0],
            [0.0, 0.0, 0.0],
            true,
        )
    end

    @testset "frechet_dist_upper_bound is symmetric and an upper bound" begin
        a = mk([(0.0, 0.0), (0.5, 0.0)])
        b = mk([(0.0, 1.0), (1.0, 1.0), (2.0, 1.0)])
        @test FD.frechet_dist_upper_bound(a, b) == FD.frechet_dist_upper_bound(b, a)
        @test FD.frechet_dist_upper_bound(a, b) >=
              FD.frechet_c_compute(a, b).leash * (1 - 1e-9)
    end
end
