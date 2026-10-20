using FrechetDist
using BasicCompGeometry
using Random
using Test

const FD = FrechetDist

include(joinpath(@__DIR__, "agreement.jl"))

@testset "Type aliases are constants" begin
    for name in
        (:DictVERType, :HeapVerticesType, :DictVerticesType, :DictHandledType, :Morphing2F)
        @test isconst(FD, name)
    end
end

@testset "Morphing validation looks at every coordinate" begin
    P = PntSeq([point(0.0, 0.0, 0.0), point(1.0, 1.0, 1.0), point(2.0, 0.0, 1.0)])
    Q = PntSeq([point(0.0, 1.0, 0.0), point(2.0, 1.0, 1.0)])
    m = FD.frechet_ve_r_compute(P, Q)
    FD.Morphing_verify_valid(m)
    bad = deepcopy(m)
    p = bad.pes[2].p
    bad.pes[2] = FD.EventPoint(
        typeof(p)(p[1], p[2], NaN),
        bad.pes[2].i,
        bad.pes[2].type,
        bad.pes[2].t,
        bad.pes[2].penalty,
    )
    @test_throws ArgumentError FD.Morphing_verify_valid(bad)
end

@testset "Polygon hierarchy" begin
    rng = MersenneTwister(11)
    P = mkseq(gen_walk(rng, 200))
    # The lock is optional (a keyword); the positional form still works.
    h = FD.ph_init(P; lock = ReentrantLock())
    h2 = FD.ph_init(P, ReentrantLock())
    h3 = FD.ph_init(P)
    @test length(h.polys) == length(h2.polys) == length(h3.polys) == 1
    FD.ph_compute_hierarchy(h)
    # The returned curve is a copy: changing it must not corrupt the hierarchy.
    a, wa = FD.ph_approx(h, 0.5)
    original = copy(a.pnts)
    a.pnts[1] = point(1.0e6, 1.0e6)
    b, wb = FD.ph_approx(h, 0.5)
    @test b.pnts == original && wb == wa
    for level in h.polys
        @test level.pnts[1] != point(1.0e6, 1.0e6)
    end
end

@testset "Size guard of the quadratic-space routines" begin
    @test FD.check_quadratic_size(100, 100, 8.2, "test") === nothing
    err = try
        FD.check_quadratic_size(10^7, 10^7, 8.2, "frechet_d_compute")
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("GiB", err.msg) && occursin("frechet_d_compute", err.msg)
end
