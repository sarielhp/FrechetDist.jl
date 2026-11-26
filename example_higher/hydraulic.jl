#!/usr/bin/env julia
# Pairwise continuous Fréchet distances of six curves in ℝ^6000: the experiment of
# Meintrup, Munteanu and Rohde ("Random projections and sampling algorithms for
# clustering of high-dimensional polygonal curves", NeurIPS 2019, Section 5) on the UCI
# hydraulic test rig data. Curve k is pressure sensor PSk, with one vertex per machine
# cycle (2205) whose coordinates are the cycle's 6000 samples. The paper reports "roughly
# three hours" for its sequential Alt-Godau implementation, without saying whether per
# distance or for all of them.
#
# Only the first CYCLES vertices and the first DIMS coordinates (default 6000) are used,
# so a run takes minutes, not hours. CYCLES may list several sizes (default 50,100):
# the code is compiled once, and the growth of the running time between the last two
# sizes gives the estimate for the full curves. If fred_hydraulic.py was run with the
# same CYCLES and DIMS, its distances and times are listed next to ours.
#
# Usage (from the repository root, after `ruby example_higher/download_hydraulic.rb`):
#   julia example_higher/hydraulic.jl [CYCLES[,CYCLES...] [DIMS]]
# Points are `Vector`s of run-time length (`GPntSeq{AnyDim}`), so the code is compiled
# once for every dimension, in seconds. (With static points, `Point{6000,Float64}`, it
# was compiled for the dimension: 11-14 min for 6000.) Times below exclude compilation.
using QuickEnv # local
using FrechetDist
using BasicCompGeometry: AnyDim, GPntSeq
using Printf

const DATA = joinpath(@__DIR__, "data", "hydraulic")
const OUTPUT = joinpath(@__DIR__, "output")
const SENSORS = 1:6
const FULL_CYCLES = 2205
# Fred stops its search within 1e-7 percent; anything beyond that and rounding is a bug.
const AGREEMENT = 1e-8

"""
    read_sensor(k, cycles, dims) -> Matrix{Float64}

The first `cycles` rows of `PSk.txt`, truncated to `dims` values, as a `dims × cycles`
matrix: column j is vertex j of the curve.
"""
function read_sensor(k::Int, cycles::Int, dims::Int)
    path = joinpath(DATA, "PS$k.txt")
    isfile(path) ||
        error("$path is missing; run `ruby example_higher/download_hydraulic.rb` first")
    M = Matrix{Float64}(undef, dims, cycles)
    rows = 0
    for line in eachline(path)
        rows == cycles && break
        rows += 1
        fields = split(line, '\t')
        length(fields) >= dims || error("$path: row $rows has $(length(fields)) values")
        for i in 1:dims
            M[i, rows] = parse(Float64, fields[i])
        end
    end
    rows == cycles || error("$path has only $rows cycles")
    return M
end

# The columns of M as `Vector` points. (`MatVecPntSeq(M)` would also do: the algorithms
# copy it into this form on entry.)
to_curve(M::Matrix{Float64}) =
    GPntSeq{AnyDim,Float64,Vector{Float64}}([M[:, j] for j in axes(M, 2)])

# Fred's results for the same truncation, as (i, j) => (distance, seconds), if present.
function read_fred(cycles::Int, dims::Int)
    path = joinpath(OUTPUT, "hydraulic_fred_$(cycles)x$(dims).csv")
    isfile(path) || return nothing
    fred = Dict{Tuple{Int,Int},Tuple{Float64,Float64}}()
    for line in Iterators.drop(eachline(path), 1)
        i, j, d, s = split(line, ',')
        fred[(parse(Int, i), parse(Int, j))] = (parse(Float64, d), parse(Float64, s))
    end
    return fred
end

# Computes one tiny instance, so the code is compiled before timing.
function warm_up(mats)
    tiny = [to_curve(M[:, 1:3]) for M in mats[1:2]]
    return frechet_c_compute(tiny[1], tiny[2]).leash
end

function pairwise(mats)
    curves = [to_curve(M) for M in mats]
    results = Tuple{Int,Int,Float64,Float64}[]
    for i in SENSORS, j in SENSORS
        i < j || continue
        seconds = @elapsed m = frechet_c_compute(curves[i], curves[j])
        push!(results, (i, j, m.leash, seconds))
    end
    return results
end

function report(results, fred)
    println(
        "pair       Fréchet (ours)    time" *
        (fred === nothing ? "" : "     Fred           time    rel. diff"),
    )
    worst = fred === nothing ? NaN : 0.0
    for (i, j, d, s) in results
        @printf("PS%d-PS%d  %16.6f %8.3f s", i, j, d, s)
        if fred !== nothing
            fd, fs = fred[(i, j)]
            rel = abs(fd - d) / d
            worst = max(worst, rel)
            @printf("  %16.6f %8.3f s  %9.2e", fd, fs, rel)
        end
        println()
    end
    total = sum(r[4] for r in results)
    @printf("total %.2f s (ours)", total)
    fred === nothing || @printf(", %.2f s (Fred)", sum(last, values(fred)))
    println()
    return total, worst
end

function write_csv(results, cycles, dims)
    mkpath(OUTPUT)
    path = joinpath(OUTPUT, "hydraulic_julia_$(cycles)x$(dims).csv")
    open(path, "w") do io
        println(io, "i,j,distance,seconds")
        for (i, j, d, s) in results
            println(io, "$i,$j,$(repr(d)),$s")
        end
    end
    println("written to $path")
end

# Pairwise distances of the first n vertices of each curve.
function run_size(mats, n::Int, dims::Int)
    results = pairwise([M[:, 1:n] for M in mats])
    println("\n$n vertices per curve")
    total, worst = report(results, read_fred(n, dims))
    write_csv(results, n, dims)
    isnan(worst) && return total, true
    ok = worst <= AGREEMENT
    @printf(
        "%s largest relative difference to Fred: %.2e\n",
        ok ? "ok    " : "FAILED",
        worst
    )
    return total, ok
end

function estimate(sizes, totals)
    if length(sizes) < 2
        @printf(
            "full curves: roughly %.0f s if time grows quadratically (not measured; pass several CYCLES)\n",
            totals[end] * (FULL_CYCLES / sizes[end])^2
        )
        return
    end
    a, b = sizes[end-1], sizes[end]
    rate = log(totals[end] / totals[end-1]) / log(b / a)
    @printf(
        "time grew as n^%.2f from %d to %d vertices; at that rate the full curves take about %.0f s\n",
        rate,
        a,
        b,
        totals[end] * (FULL_CYCLES / b)^rate
    )
end

function main(args)
    sizes = isempty(args) ? [50, 100] : parse.(Int, split(args[1], ','))
    dims = length(args) >= 2 ? parse(Int, args[2]) : 6000
    issorted(sizes) && 3 <= sizes[1] && sizes[end] <= FULL_CYCLES ||
        error("CYCLES must increase, between 3 and $FULL_CYCLES")
    println("6 curves of up to $(sizes[end]) of $FULL_CYCLES vertices, in ℝ^$dims")
    load = @elapsed mats = [read_sensor(k, sizes[end], dims) for k in SENSORS]
    println("compiling ...")
    # invokelatest (also for run_size): otherwise the code is compiled with main, before
    # the clock starts.
    compile = @elapsed Base.invokelatest(warm_up, mats)
    @printf("read in %.1f s, compiled in %.1f s\n", load, compile)
    runs = [Base.invokelatest(run_size, mats, n, dims) for n in sizes]
    println()
    estimate(sizes, first.(runs))
    return all(last, runs)
end

main(ARGS) || exit(1)
