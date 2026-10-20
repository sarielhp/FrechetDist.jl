#! /usr/bin/env julia

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

    output_polygons_to_file(list, "curves.pdf", true)
    println("Generated curves.pdf")
    output_polygons_to_file(list, "curves.png", false)
    println("Generated curves.png")
end

####################################################################

if ! isdir("output")
    mkdir("output")
end

num_args = length(ARGS);


if num_args == 0
    println("show_polygon [txt/plt file] ... [txt/plt file]")
    exit(-1)
end

plt_show(ARGS);

############################################################################
############################################################################
############################################################################
