using FrechetDist
using BasicCompGeometry
using Test

const FD = FrechetDist

curve(coords...) = PntSeq([point(c...) for c in coords])

@testset "Library code throws instead of exiting the process" begin
    # Exported validators and combinators used to call exit(-1), killing the caller.
    P = curve((0.0, 0.0), (1.0, 1.0), (2.0, 0.0))
    Q = curve((0.0, 1.0), (2.0, 1.0))
    m = FD.frechet_c_compute(P, Q)
    bad = deepcopy(m)
    bad.pes[2] = FD.EventPoint(bad.pes[2].p, bad.pes[2].i, bad.pes[2].type, 5.0, 0.0)
    @test_throws ArgumentError FD.Morphing_verify_valid(bad)

    # Chained VE morphings are not monotone; combining them used to exit(-1).
    A = curve(
        (0.26095696428881987, 0.15647047071296805),
        (0.060781372107448384, 0.5276304130131035),
        (0.5012469795915337, 0.29866073915823943),
    )
    B = curve(
        (0.042059039170194046, -0.11416283315021847),
        (-0.08663498918397372, 0.23678448083536696),
        (0.38886573287276116, -0.23119786934103614),
    )
    C = curve(
        (0.34811794631070225, 0.08228364639547436),
        (0.4837363274480386, -0.4108274804686478),
        (0.2079145053165521, -0.13887268438097466),
    )
    u = FD.frechet_ve_r_compute(A, B)
    v = FD.frechet_ve_r_compute(B, C)
    @test_throws ArgumentError FD.Morphing_combine(u, v)
end
