module TestMorphingEvents

using FrechetDist
using BasicCompGeometry
using Random
using Test

const FD = FrechetDist

function curve(rng, n)
    P = PntSeq{2,Float64}()
    x = 0.0
    for _ in 1:n
        x += rand(rng) - 0.3
        push!(P, point(x, randn(rng)))
    end
    return P
end

# Vertices passed between two consecutive events of a curve without an event at them:
# the morphing would cut the corner there.
skipped(e, f) = (e.i+1):(f.type==FD.PT_VERTEX ? f.i-1 : f.i)
corner_cuts(evs) = sum(k -> length(skipped(evs[k], evs[k+1])), 1:(length(evs)-1); init = 0)

@testset "Morphings visit every vertex" begin
    # Events computed a few ulps off a vertex used to be recorded inside the next edge
    # (or at t = 1 of the previous one), so the vertex between two events was skipped.
    fns = [
        ("frechet_c_compute", (P, Q) -> frechet_c_compute(P, Q)),
        ("frechet_c_approx", (P, Q) -> frechet_c_approx(P, Q, 1.1)),
        (
            "frechet_mono_via_refinement",
            (P, Q) -> frechet_mono_via_refinement(P, Q, 1.001)[1],
        ),
        ("frechet_ve_r_compute", (P, Q) -> frechet_ve_r_compute(P, Q)),
    ]
    for (name, f) in fns
        @testset "$name" begin
            rng = Xoshiro(3)
            cuts = 0
            for _ in 1:40
                m = f(curve(rng, 30), curve(rng, 30))
                cuts += corner_cuts(m.pes) + corner_cuts(m.qes)
            end
            @test cuts == 0
        end
    end
end

end # module
