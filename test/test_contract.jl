module TestContract

using FrechetDist
using BasicCompGeometry
using Test

const FD = FrechetDist

mk(coords) = PntSeq([point(c...) for c in coords])
const two = mk([(0.0, 0.0), (1.0, 0.0)])
const three = mk([(0.0, 1.0), (1.0, 1.5), (2.0, 1.0)])
const one = mk([(0.0, 0.0)])
const empty = PntSeq2F()

# name => (function of two curves, minimum number of vertices)
const entry_points = Dict(
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

    @testset "Element types of single-curve entry points" begin
        # Only Float64 coordinates are supported; others get an ArgumentError, not a
        # MethodError or InexactError from deep inside.
        coords = ((0, 0), (1, 1), (2, 0), (3, 1))
        for T in (Float32, Int)
            P = PntSeq([Point{2,T}(c...) for c in coords])
            calls = (
                () -> FD.frechet_width_approx(P),
                () -> FD.frechet_offsets(P, [1, 2, 4]),
                () -> FD.frechet_palette(P),
                () -> FD.frechet_palette_approx(P, T[0, 0, 0, 0], oneunit(T)),
                () -> FD.frechet_simplify_to_width(P, oneunit(T)),
                () -> FD.frechet_simplify_w_exp(P, oneunit(T)),
                () -> FD.frechet_c_mono_approx_subcurve(P, P, collect(1:4)),
                () -> FD.frechet_ve_r_compute_range(P, P, oneunit(T)),
                () -> FD.FEVER_compute_range(P, P, oneunit(T)),
                () -> FD.Morphing_empty(P, P),
            )
            for f in calls
                @test outcome((_, _) -> f(), nothing, nothing) == :ArgumentError
            end
        end
    end

    @testset "Widths and bounds accept any real number" begin
        P = mk([(0.0, 0.0), (1.0, 0.3), (2.0, 0.0), (3.0, 0.4), (4.0, 0.0)])
        # PntSeq has no ==; compare the points and the remaining results.
        same(a, b) = a[1].pnts == b[1].pnts && a[2:end] == b[2:end]
        @test same(FD.frechet_simplify_to_width(P, 1), FD.frechet_simplify_to_width(P, 1.0))
        @test same(FD.frechet_simplify_w_exp(P, 1), FD.frechet_simplify_w_exp(P, 1.0))
        @test FD.frechet_ve_r_compute_range(P, three, 10) ==
              FD.frechet_ve_r_compute_range(P, three, 10.0)
        @test FD.FEVER_compute_range(P, three, 10) == FD.FEVER_compute_range(P, three, 10.0)
        plt = FD.frechet_palette(P)
        @test same(
            FD.frechet_palette_approx(P, plt, 1),
            FD.frechet_palette_approx(P, plt, 1.0),
        )
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

@testset "Sampling and subcurve functions" begin
    # Zero-length curves used to give InexactError (ceil(NaN)); n < 1 silently returned a value.
    point_curve = mk([(1.0, 1.0), (1.0, 1.0), (1.0, 1.0)])
    other_point = mk([(1.0, 3.0), (1.0, 3.0)])
    for lopt in (true, false)
        @test FD.frechet_d_r_compute_sample(point_curve, other_point, 5, lopt).leash ≈ 2.0
        @test FD.frechet_d_r_compute_sample(two, three, 10, lopt).leash > 0
        @test_throws ArgumentError FD.frechet_d_r_compute_sample(two, three, 0, lopt)
        @test_throws ArgumentError FD.frechet_d_r_compute_sample(two, three, -3, lopt)
    end

    # Vertex index lists must be in range and non-decreasing.
    curve = mk([(0.0, 0.0), (1.0, 0.5), (2.0, 0.0), (3.0, 0.5), (4.0, 0.0)])
    @test length(FD.frechet_offsets(curve, [1, 3, 5])) == 2
    @test_throws ArgumentError FD.frechet_offsets(curve, [1, 6])
    @test_throws ArgumentError FD.frechet_offsets(curve, [0, 3])
    @test_throws ArgumentError FD.frechet_offsets(curve, [5, 3, 1])
    simplified = mk([(0.0, 0.0), (2.0, 0.0), (4.0, 0.0)])
    @test FD.frechet_c_mono_approx_subcurve(curve, simplified, [1, 3, 5])[1].leash >= 0
    @test_throws ArgumentError FD.frechet_c_mono_approx_subcurve(curve, simplified, [1, 3])
    @test_throws ArgumentError FD.frechet_c_mono_approx_subcurve(
        curve,
        simplified,
        [1, 3, 9],
    )
    @test_throws ArgumentError FD.frechet_c_mono_approx_subcurve(
        curve,
        simplified,
        [1, 5, 3],
    )
end


@testset "Refinement loops give up instead of hanging" begin
    saved = FD.MAX_REFINEMENT_ROUNDS[]
    try
        FD.MAX_REFINEMENT_ROUNDS[] = 0
        @test_throws ErrorException FD.frechet_c_compute(two, three)
        @test_throws ErrorException FD.frechet_c_approx(two, three, 1.5)
        @test_throws ErrorException FD.frechet_mono_via_refinement(two, three)
    finally
        FD.MAX_REFINEMENT_ROUNDS[] = saved
    end
    @test FD.frechet_c_compute(two, three).leash > 0
end

end # module
