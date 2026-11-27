#!/usr/bin/env julia
# Exact continuous Fréchet distances (frechet_c_compute_exact) between all pairs of a
# set of curves, and how far the other results are from them: those of
# frechet_c_compute (curves.jl, hydraulic.jl) and of Fred (fred_curves.py,
# fred_hydraulic.py), if their output files are present.
#
# Usage: julia [-t THREADS] example_higher/exact_curves.jl SET [FIRST LAST]
#   SET = pems or kth (curves FIRST to LAST of data/curves/SET, default 1 to 14), or
#   hydraulic (the six full curves PS1..PS6, 2205 vertices in ℝ^6000).
# Writes output/TAG_exact.csv (a,b,distance,seconds); TAG is e.g. pems_1-14 or
# hydraulic_2205x6000, as for the other scripts. The tables of frechet_c_compute_exact
# are built on all of Julia's threads (the result does not depend on their number).
using QuickEnv # local
using FrechetDist
using BasicCompGeometry: AnyDim, GPntSeq
using Printf

const CURVES = joinpath(@__DIR__, "data", "curves")
const HYDRAULIC = joinpath(@__DIR__, "data", "hydraulic")
const OUTPUT = joinpath(@__DIR__, "output")
const Curve = GPntSeq{AnyDim,Float64,Vector{Float64}}
# The results compared are those at the tolerance FRECHET_TOL (as curves.jl,
# hydraulic.jl and fred_curves.py name their files; default 1e-9).
const TOL_ENV = get(ENV, "FRECHET_TOL", "1e-9")
const TOL_SUFFIX = parse(Float64, TOL_ENV) == 1e-9 ? "" : "_tol" * TOL_ENV

# A prepared set (prepare_pems.rb, prepare_kth.rb): rows FIRST..LAST of index.csv, each
# a file of dims x vertices Float64 values.
function read_set(set::String, first::Int, last::Int)
    rows = [
        split(l, ',') for
        l in Iterators.drop(eachline(joinpath(CURVES, set, "index.csv")), 1)
    ]
    names, curves = String[], Curve[]
    for r in rows[first:last]
        M = Matrix{Float64}(undef, parse(Int, r[3]), parse(Int, r[4]))
        read!(joinpath(CURVES, set, r[5]), M)
        push!(names, r[1])
        push!(curves, Curve([M[:, j] for j in axes(M, 2)]))
    end
    return names, curves
end

# The hydraulic curves: row j of PSk.txt (6000 tab-separated values) is vertex j.
function read_hydraulic()
    curves = map(1:6) do k
        path = joinpath(HYDRAULIC, "PS$k.txt")
        isfile(path) ||
            error("$path is missing; run example_higher/download_hydraulic.rb")
        Curve([parse.(Float64, split(line, '\t')) for line in eachline(path)])
    end
    return string.(1:6), curves
end

# Another method's results for the same pairs, as (a, b) => (distance, seconds).
function read_results(path::String)
    isfile(path) || return nothing
    out = Dict{Tuple{String,String},Tuple{Float64,Float64}}()
    for line in Iterators.drop(eachline(path), 1)
        a, b, d, s = split(line, ',')
        out[(a, b)] = (parse(Float64, d), parse(Float64, s))
    end
    return out
end

# The results file of another method: SET_FIRST-LAST_METHOD.csv, or for the hydraulic
# curves hydraulic_METHOD_2205x6000.csv (the names hydraulic.jl and fred_hydraulic.py use).
function other_file(tag::String, method::String)
    set, rest = split(tag, "_"; limit = 2)
    name =
        set == "hydraulic" ? "$(set)_$(method)_$rest$TOL_SUFFIX" :
        "$(tag)$(TOL_SUFFIX)_$method"
    return joinpath(OUTPUT, "$name.csv")
end

function pairwise(names, curves)
    # Compiled before the clock starts.
    frechet_c_compute_exact(Curve(curves[1].pnts[1:3]), Curve(curves[2].pnts[1:3]))
    results = Tuple{String,String,Float64,Float64}[]
    for i in eachindex(curves), j in (i+1):length(curves)
        seconds = @elapsed d = frechet_c_compute_exact(curves[i], curves[j])
        push!(results, (names[i], names[j], d, seconds))
        @printf("%-14s %-14s %22.16g %9.3f s\n", names[i], names[j], d, seconds)
        flush(stdout)
    end
    return results
end

# How far another method's distances are from the exact ones.
function compare(label::String, other, results)
    other === nothing && return println("$label: no results file")
    rel = [(other[(a, b)][1] - d) / d for (a, b, d, _) in results]
    total = sum(other[(a, b)][2] for (a, b, _, _) in results)
    worst = argmax(abs, rel)
    @printf(
        "%-18s %9.1f s; worst relative error %+.2e; off by more than 1e-8: %d (%d below)\n",
        label,
        total,
        worst,
        count(>(1e-8) ∘ abs, rel),
        count(<(-1e-8), rel)
    )
end

function main(args)
    set = isempty(args) ? "pems" : args[1]
    if set == "hydraulic"
        tag = "hydraulic_2205x6000"
        names, curves = read_hydraulic()
    else
        first, last = length(args) >= 3 ? parse.(Int, args[2:3]) : (1, 14)
        tag = "$(set)_$first-$last"
        names, curves = read_set(set, first, last)
    end
    println(
        "$tag: $(length(curves)) curves in ℝ^$(length(curves[1].pnts[1])), ",
        "$(Threads.nthreads()) threads",
    )
    results = Base.invokelatest(pairwise, names, curves)
    mkpath(OUTPUT)
    open(joinpath(OUTPUT, "$(tag)_exact.csv"), "w") do io
        println(io, "a,b,distance,seconds")
        foreach(r -> println(io, join(r, ',')), results)
    end
    @printf(
        "\n%-18s %9.1f s for %d pairs\n",
        "exact",
        sum(r[4] for r in results),
        length(results)
    )
    compare("frechet_c_compute", read_results(other_file(tag, "julia")), results)
    compare("Fred", read_results(other_file(tag, "fred")), results)
end

main(ARGS)
