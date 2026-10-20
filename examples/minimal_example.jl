# julia

push!(LOAD_PATH, pwd()*"/src/")

using FrechetDist
using BasicCompGeometry

include("graphics.jl")


function min_example()
    P = PntSeq2F()
    Q = PntSeq2F()

    push!(P, point(0.0, 1.0), point(10.0, 1.0), point(10.0, 2.0), point(20.0, 2.0))
    push!(Q, point(0.0, -1.0), point(10.0, -2.0), point(20.0, 3.0))

    return frechet_c_compute(P, Q)
end


mw = min_example()
output_morphing(mw, "morphing.pdf");
