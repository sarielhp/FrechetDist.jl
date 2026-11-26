#!/usr/bin/env julia
# Static points (`SVector{D}`) against run-time-length points (`Vector{Float64}`) for the
# kernels the Fréchet code spends its time in, across dimensions: which point type is
# faster, and where the crossover lies.
#
# Static points carry the dimension in their type. The compiler unrolls every loop over
# the coordinates and keeps points out of the heap, which is ideal in the plane; but the
# code is compiled again for each dimension, the compile time grows superlinearly with
# it, and the unrolled code runs slower than a plain loop once D is large. Vector points
# compile once for all dimensions and run a vectorized loop, at the price of a heap
# allocation for every point they create.
#
# Kernels (FrechetDist's `fp_` functions for static points; the loop versions below,
# with the same arithmetic in the same order, for vector points):
#   dist      distance between two points
#   dps       distance from a point to a segment (fused into one pass for vectors)
#   comb      a convex combination of two points, which creates a point, and its
#             distance to a third
# Each kernel runs on every vertex/vertex or vertex/edge pair of two random walks of
# N = 100 vertices. Each (point type, D) runs in its own process, so the reported
# compile time is that of the first call in a fresh session.
#
# Result (AMD Ryzen 7 5700G, Julia 1.13.3; the README has the table): static points are
# faster up to D = 16, it is close from 32 to 96, and vector points are faster from
# D = 128 on (distances 2-3 times, point-to-segment 1.7-4.5 times), except for creating
# points, where the allocation keeps static points ahead up to D = 512. Static kernels
# compile in 0.14 s at small D, 4.5 s at 1024, 13.5 s at 2048; vector ones in 0.2 s.
#
# Usage (from the repository root):
#   julia example_higher/point_types.jl [D,D,...]
# Prints a table and writes output/point_types.csv.
using QuickEnv # local
using FrechetDist
using BasicCompGeometry: Point
using Printf
using Random

const FD = FrechetDist
const N = 100
const KERNELS = ("dist", "dps", "comb")
const DEFAULT_DIMS =
    [2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, 512, 1024, 2048]
const OUTPUT = joinpath(@__DIR__, "output")

# --- vector-point kernels: the arithmetic of the fp_ versions, without temporary vectors

@inline function vec_dist_sq(p::Vector{Float64}, q::Vector{Float64})
    s = 0.0
    @inbounds for i in eachindex(p, q)
        @fastmath s += (p[i] - q[i])^2 # julia-audit: allow-fastmath (as fp_dist_sq)
    end
    return s
end

@inline vec_dist(p::Vector{Float64}, q::Vector{Float64}) = sqrt(vec_dist_sq(p, q))

@inline function vec_convex_comb(p::Vector{Float64}, q::Vector{Float64}, t::Float64)
    t == 0.0 && return p  # exact, as fp_convex_comb
    t == 1.0 && return q  # exact, as above
    o = similar(p)
    @inbounds for i in eachindex(p, q)
        @fastmath o[i] = p[i] + (q[i] - p[i]) * t # julia-audit: allow-fastmath (as fp_convex_comb)
    end
    return o
end

# Squared distance from qr to the point s_p + t (s_q - s_p): fp_dist_point_segment's
# fallback when a t^2 + b t + c cancels.
function vec_dist_sq_at(qr, s_p, s_q, t::Float64)
    val = 0.0
    @inbounds for i in eachindex(qr, s_p, s_q)
        r = (s_p[i] - qr[i]) + t * (s_q[i] - s_p[i])
        val += r * r
    end
    return val
end

# fp_dist_point_segment, with <s_p - qr, s_q - s_p> and |s_p - qr|^2 in one pass.
function vec_dist_point_segment(
    qr::Vector{Float64},
    s_p::Vector{Float64},
    s_q::Vector{Float64},
)
    a = vec_dist_sq(s_p, s_q)
    bh = 0.0
    c = 0.0
    @inbounds for i in eachindex(qr, s_p, s_q)
        d = s_p[i] - qr[i]
        bh += d * (s_q[i] - s_p[i])
        c += d * d
    end
    b = 2.0 * bh
    val = c
    if a != 0.0
        t = -b / (2.0 * a)
        if t > 1
            val = vec_dist_sq(qr, s_q)
        elseif t >= 0
            val = a * t^2 + b * t + c
            val < 1.0e-8 * c && (val = vec_dist_sq_at(qr, s_p, s_q, t))
        end
    end
    val <= FD.FP_ZERO_SQ * max(a, c) && return 0.0
    return sqrt(val)
end

# --- workloads: all pairs of two curves; the sum keeps the work from being optimized away

function work_dist(P, Q, dist)
    s = 0.0
    for p in P, q in Q
        s += dist(p, q)
    end
    return s
end

function work_dps(P, Q, dps)
    s = 0.0
    for p in P, j in 1:(length(Q)-1)
        s += dps(p, Q[j], Q[j+1])
    end
    return s
end

# The new point is used in full (a distance), or the compiler computes only what is read.
function work_comb(P, Q, comb, dist)
    s = 0.0
    for i in eachindex(P), j in 1:(length(Q)-1)
        s += dist(comb(Q[j], Q[j+1], i / (length(P) + 1)), P[i])
    end
    return s
end

function random_walk(D::Int, seed::Int)
    rng = MersenneTwister(seed)
    cur = zeros(D)
    walk = Vector{Vector{Float64}}()
    for _ in 1:N
        cur = cur .+ (rand(rng, D) .- 0.5)
        push!(walk, copy(cur))
    end
    return walk
end

# (seconds spent running, seconds per run, value): the best of five batches, each long
# enough (50 ms) to time reliably. Compiling `measure` for `f` compiles `f` as well,
# before any timer here starts; `timed` measures that from the outside.
function measure(f)
    start = time_ns()
    v = f()
    once = @elapsed f()
    reps = max(1, ceil(Int, 0.05 / max(once, 1.0e-9)))
    best = Inf
    for _ in 1:5
        best = min(best, @elapsed(for _ in 1:reps
            v = f()
        end) / reps)
    end
    return (time_ns() - start) / 1.0e9, best, v
end

# (compile seconds, seconds per run, value). `invokelatest` keeps the compilation of
# `measure(f)` inside the timed call.
function timed(f)
    total = @elapsed (running, run, v) = Base.invokelatest(measure, f)
    return max(total - running, 0.0), run, v
end

function kernel_runs(P, Q, static::Bool)
    if static
        return (
            () -> work_dist(P, Q, FD.fp_dist),
            () -> work_dps(P, Q, FD.fp_dist_point_segment),
            () -> work_comb(P, Q, FD.fp_convex_comb, FD.fp_dist),
        )
    end
    return (
        () -> work_dist(P, Q, vec_dist),
        () -> work_dps(P, Q, vec_dist_point_segment),
        () -> work_comb(P, Q, vec_convex_comb, vec_dist),
    )
end

static_points(walk, ::Val{D}) where {D} = [Point{D,Float64}(w) for w in walk]

# One (point type, D) measurement; prints "kernel compile_s ns_per_call value" lines.
function child(kind::String, D::Int)
    P, Q = random_walk(D, 1), random_walk(D, 2)
    static = kind == "static"
    if static
        P, Q = static_points(P, Val(D)), static_points(Q, Val(D))
    end
    for (name, f) in zip(KERNELS, kernel_runs(P, Q, static))
        calls = name == "dist" ? N * N : N * (N - 1)
        compile, run, v = timed(f)
        println(name, " ", compile, " ", 1.0e9 * run / calls, " ", repr(v))
    end
end

function run_child(kind::String, D::Int)
    cmd = `$(Base.julia_cmd()) --project=$(@__DIR__) $(@__FILE__) --child $kind $D`
    rows = Dict{String,Tuple{Float64,Float64,String}}()
    for line in eachline(cmd)
        name, compile, ns, v = split(line)
        rows[name] = (parse(Float64, compile), parse(Float64, ns), v)
    end
    return rows
end

function print_row(D, s, v)
    @printf("%6d", D)
    for k in KERNELS
        @printf(
            "  %9.1f %9.1f %5.2f%s",
            s[k][2],
            v[k][2],
            s[k][2] / v[k][2],
            s[k][3] == v[k][3] ? " " : "*"
        )
    end
    @printf("  %8.2f %6.2f\n", sum(s[k][1] for k in KERNELS), sum(v[k][1] for k in KERNELS))
end

function write_csv(results)
    mkpath(OUTPUT)
    path = joinpath(OUTPUT, "point_types.csv")
    open(path, "w") do io
        println(io, "D,kernel,points,compile_s,ns_per_call,value")
        for (D, s, v) in results, k in KERNELS, (kind, r) in (("static", s), ("vector", v))
            println(io, "$D,$k,$kind,$(r[k][1]),$(r[k][2]),$(r[k][3])")
        end
    end
    println("written to $path")
end

# The smallest tested D from which on vector points are faster for the kernel.
function crossover(results, k)
    i = findlast(r -> r[2][k][2] <= r[3][k][2], results)
    i === nothing && return results[1][1]
    i == length(results) && return nothing
    return results[i+1][1]
end

function main(args)
    dims = isempty(args) ? DEFAULT_DIMS : parse.(Int, split(args[1], ','))
    println("ns per call, static / vector points, and the ratio (> 1: vectors are faster);")
    println("* marks a kernel whose sums differ. Last columns: compile s, static / vector.")
    @printf("%6s  %-27s  %-27s  %-27s  %s\n", "D", "dist", "dps", "comb", "compile")
    results = []
    for D in dims
        s, v = run_child("static", D), run_child("vector", D)
        print_row(D, s, v)
        push!(results, (D, s, v))
    end
    for k in KERNELS
        c = crossover(results, k)
        println("$k: vectors faster from D = ", c === nothing ? "(not in range)" : c)
    end
    write_csv(results)
end

if length(ARGS) == 3 && ARGS[1] == "--child"
    child(ARGS[2], parse(Int, ARGS[3]))
else
    main(ARGS)
end
