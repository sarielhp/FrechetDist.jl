#! julia

push!(LOAD_PATH, pwd()*"/src/")

using Cairo
using FrechetDist
using BasicCompGeometry
import BasicCompGeometry: translate!, bound!
import FrechetDist: read_file
using Parameters

include("graphics.jl")

function plt_show(ARGS)
    num_args = length(ARGS)

    if num_args == 0
        return
    end


    list = Vector{PntSeq2F}()
    bb = BBox2F()
    for i = 1:num_args
        println("Reading: ", ARGS[i])
        poly_a = read_file(ARGS[i])
        push!(list, poly_a)
        bound!(bb, poly_a)
    end
    p = bottom_left(bb)
    for poly in list
        translate!(poly, -p)
    end

    out = Vector{PntSeq2F}()
    P = list[1]
    push!(out, P)

    r = dist(first(P), last(P)) / 30.0
    Q = hausdorff_simplify(P, r)[1]
    push!(out, Q)

    output_polygons_to_file(out, "h.pdf", true)
    output_polygons_to_file(out, "h.png", false)
    println("Generated h.png h.pdf")
end

####################################################################

if ! isdir("output")
    mkdir("output")
end

num_args = length(ARGS);

if num_args == 0
    println("simplify_hausdorff [txt/plt file] ... [txt/plt file]")
    exit(-1)
end

plt_show(ARGS);

############################################################################
############################################################################
############################################################################
