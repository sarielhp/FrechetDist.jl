using FrechetDist
using BasicCompGeometry
using Random
using Test

const FD = FrechetDist

include(joinpath(@__DIR__, "agreement.jl"))

# Everything the code writes to stdout while running `f`.
function stdout_of(f)
    return mktemp() do path, io
        redirect_stdout(io) do
            try
                f()
            catch
            end
        end
        flush(io)
        read(path, String)
    end
end

@testset "Library code does not print" begin
    rng = MersenneTwister(3)
    P = mkseq(gen_walk(rng, 40))
    Q = mkseq(gen_walk(rng, 35; step = 2.0))
    m = FD.frechet_c_compute(P, Q)
    plt = FD.frechet_palette(P)
    calls = Dict(
        "frechet_palette_level" => () -> FD.frechet_palette_level(P, plt, 2),
        "frechet_palette_approx" => () -> FD.frechet_palette_approx(P, plt, 1.0),
        "simplify_morphing_sensitive" => () -> FD.simplify_morphing_sensitive(m, 4.0),
        "frechet_simplify_to_cardin" => () -> FD.frechet_simplify_to_cardin(P, 8),
        "frechet_c_compute" => () -> FD.frechet_c_compute(P, Q),
        "frechet_c_approx" => () -> FD.frechet_c_approx(P, Q, 1.1),
        "frechet_ve_r_compute_range" => () -> FD.frechet_ve_r_compute_range(P, Q, 0.01),
        "Morphing_combine (non-monotone)" =>
            () -> FD.Morphing_combine(
                FD.frechet_ve_r_compute(P, Q),
                FD.frechet_ve_r_compute(Q, P),
            ),
        "invalid input" => () -> FD.frechet_c_compute(PntSeq2F(), Q),
        "NaN input" =>
            () -> FD.frechet_ve_r_compute(mkseq([[0.0, 0.0], [NaN, 1.0]]), Q),
    )
    for (name, f) in sort(collect(calls); by = first)
        @test (name, stdout_of(f)) == (name, "")
    end
end
