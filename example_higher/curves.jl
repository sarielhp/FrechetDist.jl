#!/usr/bin/env julia
# Pairwise continuous Fréchet distances among curves of a prepared set: curves FIRST to
# LAST (rows of data/curves/SET/index.csv, written by prepare_pems.rb or
# prepare_kth.rb), with `Vector` points (GPntSeq{AnyDim}), compiled once for every
# dimension. If fred_curves.py was run on the same selection, Fred's distances and times
# are listed next to ours, with the pairs where the two differ by more than 1e-8 relative
# (Fred is sometimes too low; check_frechet.jl decides independently which is right).
# Also reports how often a curve's nearest neighbour has the same label.
#
# Usage: julia example_higher/curves.jl SET [FIRST [LAST]]     (default curves 1 to 14)
# The environment variable FRECHET_TOL sets the tolerance of frechet_c_compute (default
# here 1e-9, comparable to the 1e-7 percent Fred is run with; the package default is
# 1e-5); fred_curves.py reads the same variable.
using QuickEnv # local
using FrechetDist
using BasicCompGeometry: AnyDim, GPntSeq
using Printf

const CURVES = joinpath(@__DIR__, "data", "curves")
const OUTPUT = joinpath(@__DIR__, "output")
# Fred stops its search within 1e-7 percent, so both should agree to about 1e-8.
const AGREEMENT = 1e-8
const TOL = parse(Float64, get(ENV, "FRECHET_TOL", "1e-9"))
# Results at another tolerance than the default go to files with it in their names
# (`_tol1e-5` for FRECHET_TOL=1e-5), next to those at the default.
const TOL_SUFFIX = TOL == 1e-9 ? "" : "_tol" * ENV["FRECHET_TOL"]

const Curve = GPntSeq{AnyDim,Float64,Vector{Float64}}

"The rows of `index.csv` of a set: (name, label, dims, vertices, file)."
function read_index(set::String)
    path = joinpath(CURVES, set, "index.csv")
    isfile(path) || error("$path is missing; run the set's prepare script first")
    rows = map(Iterators.drop(eachline(path), 1)) do line
        name, label, dims, vertices, file = split(line, ',')
        (; name, label, dims = parse(Int, dims), vertices = parse(Int, vertices), file)
    end
    return rows
end

# The file holds dims x vertices Float64 values, vertex after vertex (little-endian,
# the byte order of the machines this runs on).
function read_curve(set::String, row)
    M = Matrix{Float64}(undef, row.dims, row.vertices)
    read!(joinpath(CURVES, set, row.file), M)
    return Curve([M[:, j] for j in axes(M, 2)])
end

# Fred's results for the same selection, as (a, b) => (distance, seconds), if present.
function read_fred(tag::String)
    path = joinpath(OUTPUT, "$(tag)_fred.csv")
    isfile(path) || return nothing
    fred = Dict{Tuple{String,String},Tuple{Float64,Float64}}()
    for line in Iterators.drop(eachline(path), 1)
        a, b, d, s = split(line, ',')
        fred[(a, b)] = (parse(Float64, d), parse(Float64, s))
    end
    return fred
end

warm_up(curves) =
    frechet_c_compute(Curve(curves[1].pnts[1:3]), Curve(curves[2].pnts[1:3])).leash

function pairwise(names, curves)
    results = Tuple{String,String,Float64,Float64}[]
    for i in eachindex(curves), j in eachindex(curves)
        i < j || continue
        seconds = @elapsed m = frechet_c_compute(curves[i], curves[j], true, TOL)
        push!(results, (names[i], names[j], m.leash, seconds))
    end
    return results
end

function report(results, fred)
    worst = 0.0
    differ = fred_lower = 0
    for (a, b, d, s) in results
        line = @sprintf("%-14s %-14s %14.6f %8.3f s", a, b, d, s)
        if fred !== nothing
            fd, fs = fred[(a, b)]
            rel = abs(fd - d) / max(d, eps())
            worst = max(worst, rel)
            if rel > AGREEMENT
                differ += 1
                fred_lower += fd < d
            end
            line *= @sprintf("  %14.6f %8.3f s  %9.2e", fd, fs, rel)
        end
        length(results) <= 40 && println(line)
    end
    total = sum(r[4] for r in results)
    @printf("%d pairs: %.2f s (ours)", length(results), total)
    fred === nothing || @printf(", %.2f s (Fred)", sum(last, values(fred)))
    println()
    fred === nothing && return
    @printf(
        "largest relative difference to Fred: %.2e; %d pairs differ by more than %.0e",
        worst,
        differ,
        AGREEMENT
    )
    differ > 0 && @printf(" (Fred lower in %d; see check_frechet.jl)", fred_lower)
    println()
end

# How many curves have a nearest neighbour (by Fréchet distance) with the same label.
function nearest_same_label(results, names, labels)
    dist = Dict{Tuple{String,String},Float64}()
    for (a, b, d, _) in results
        dist[(a, b)] = d
        dist[(b, a)] = d
    end
    label = Dict(zip(names, labels))
    same = count(names) do a
        nearest = argmin(b -> dist[(a, b)], filter(!=(a), names))
        label[nearest] == label[a]
    end
    @printf("nearest neighbour has the same label: %d of %d curves\n", same, length(names))
end

function write_csv(results, tag)
    mkpath(OUTPUT)
    path = joinpath(OUTPUT, "$(tag)_julia.csv")
    open(path, "w") do io
        println(io, "a,b,distance,seconds")
        foreach(r -> println(io, "$(r[1]),$(r[2]),$(repr(r[3])),$(r[4])"), results)
    end
    println("written to $path")
end

function main(args)
    isempty(args) && error("usage: julia curves.jl SET [FIRST [LAST]]")
    set = args[1]
    index = read_index(set)
    first_k = length(args) >= 2 ? parse(Int, args[2]) : 1
    last_k = length(args) >= 3 ? parse(Int, args[3]) : min(14, length(index))
    1 <= first_k < last_k <= length(index) ||
        error("need 1 <= FIRST < LAST <= $(length(index))")
    rows = index[first_k:last_k]
    load = @elapsed curves = [read_curve(set, r) for r in rows]
    @printf(
        "%s: curves %d-%d, %d to %d vertices in R^%d, read in %.1f s\n",
        set,
        first_k,
        last_k,
        minimum(length, curves),
        maximum(length, curves),
        rows[1].dims,
        load
    )
    # invokelatest: otherwise the code is compiled with main, before the clock starts.
    compile = @elapsed Base.invokelatest(warm_up, curves)
    @printf("compiled in %.1f s\n", compile)
    names = [String(r.name) for r in rows]
    results = Base.invokelatest(pairwise, names, curves)
    tag = "$(set)_$(first_k)-$(last_k)$TOL_SUFFIX"
    fred = read_fred(tag)
    report(results, fred)
    nearest_same_label(results, names, [r.label for r in rows])
    write_csv(results, tag)
end

main(ARGS)
