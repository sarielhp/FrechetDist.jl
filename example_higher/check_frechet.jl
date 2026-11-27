#!/usr/bin/env julia
# Independent check of the distances computed by curves.jl (and fred_curves.py): the
# textbook free-space decision procedure of Alt and Godau (1995), written from scratch,
# sharing no code with FrechetDist, and a binary search over it to 1 ulp. For every pair
# in output/SET_FIRST-LAST_julia.csv (and _fred.csv, if present) prints the true
# distance and the relative error of each answer, and lists the answers that are off by
# more than 1e-8 or below the true distance (an answer below it is wrong: no matching
# achieves it). O(n m d) per decision, about 60 decisions per pair.
#
# Usage: julia example_higher/check_frechet.jl SET FIRST LAST
using Printf

const CURVES = joinpath(@__DIR__, "data", "curves")
const OUTPUT = joinpath(@__DIR__, "output")
const REPORT = 1e-8
const EMPTY = (1.0, 0.0)

isempty_iv(iv) = iv[1] > iv[2]

# The t in [0, 1] with |a + t (b - a) - c| <= eps, as an interval (lo, hi).
function free_interval(c, a, b, eps)
    A = 0.0
    B = 0.0
    C = 0.0
    @inbounds for k in eachindex(c)
        d = b[k] - a[k]
        e = a[k] - c[k]
        A += d * d
        B += e * d
        C += e * e
    end
    C -= eps * eps
    A == 0.0 && return C <= 0 ? (0.0, 1.0) : EMPTY
    disc = B * B - A * C
    disc < 0 && return EMPTY
    r = sqrt(disc)
    return (max(0.0, (-B - r) / A), min(1.0, (-B + r) / A))
end

# Reachable parts of the first column of left edges (or first row of bottom edges).
function first_line!(R, p, Q, eps)
    for j in 1:(length(Q)-1)
        iv = free_interval(p, Q[j], Q[j+1], eps)
        iv[1] == 0.0 || return
        R[j] = iv
        iv[2] == 1.0 || return
    end
end

# Is the Fréchet distance of P and Q at most eps? Cell (i, j) spans P[i]P[i+1] x
# Q[j]Q[j+1]; L[i, j] is the reachable part of its left edge (P[i] against the Q
# segment), Bo[i, j] that of its bottom edge (Q[j] against the P segment).
function decide(P, Q, eps)
    n, m = length(P), length(Q)
    sum(abs2, P[1] - Q[1]) <= eps^2 && sum(abs2, P[n] - Q[m]) <= eps^2 || return false
    L = fill(EMPTY, n, m - 1)
    Bo = fill(EMPTY, n - 1, m)
    first_line!(view(L, 1, :), P[1], Q, eps)
    first_line!(view(Bo, :, 1), Q[1], P, eps)
    for i in 1:(n-1), j in 1:(m-1)
        l, b = L[i, j], Bo[i, j]
        top = free_interval(Q[j+1], P[i], P[i+1], eps)
        right = free_interval(P[i+1], Q[j], Q[j+1], eps)
        if !isempty_iv(l)
            Bo[i, j+1] = top
        elseif !isempty_iv(b)
            Bo[i, j+1] = (max(top[1], b[1]), top[2])
        end
        if !isempty_iv(b)
            L[i+1, j] = right
        elseif !isempty_iv(l)
            L[i+1, j] = (max(right[1], l[1]), right[2])
        end
    end
    a, b = L[n, m-1], Bo[n-1, m]
    return (!isempty_iv(a) && a[2] == 1.0) || (!isempty_iv(b) && b[2] == 1.0)
end

# The smallest eps (to 1 ulp) that decide accepts; `guess` is any answer to start from.
function distance(P, Q, guess)
    lo, hi = 0.0, 2 * guess + 1.0
    while true
        mid = (lo + hi) / 2
        (mid == lo || mid == hi) && return hi
        decide(P, Q, mid) ? (hi = mid) : (lo = mid)
    end
end

function read_curve(set, file, dims)
    M = Matrix{Float64}(undef, dims, filesize(joinpath(CURVES, set, file)) ÷ (8 * dims))
    read!(joinpath(CURVES, set, file), M)
    return [M[:, j] for j in axes(M, 2)]
end

# (a, b) => distance from a results file, or nothing if it does not exist.
function read_results(path)
    isfile(path) || return nothing
    return Dict(
        (String(f[1]), String(f[2])) => parse(Float64, f[3]) for
        f in (split(l, ',') for l in Iterators.drop(eachline(path), 1))
    )
end

function summarize(name, errors)
    off = count(e -> abs(e) > REPORT, errors)
    low = count(e -> e < -REPORT, errors)
    @printf(
        "%-5s worst relative error %.2e; %d of %d off by more than %.0e, %d below the true distance\n",
        name,
        maximum(abs, errors),
        off,
        length(errors),
        REPORT,
        low
    )
end

function main(args)
    length(args) == 3 || error("usage: julia check_frechet.jl SET FIRST LAST")
    set, tag = args[1], "$(args[1])_$(args[2])-$(args[3])"
    ours = read_results(joinpath(OUTPUT, "$(tag)_julia.csv"))
    ours === nothing && error("output/$(tag)_julia.csv is missing; run curves.jl first")
    fred = read_results(joinpath(OUTPUT, "$(tag)_fred.csv"))
    index = Dict(
        String(f[1]) => (parse(Int, f[3]), String(f[5])) for f in (
            split(l, ',') for
            l in Iterators.drop(eachline(joinpath(CURVES, set, "index.csv")), 1)
        )
    )
    curve(name) = read_curve(set, index[name][2], index[name][1])
    errs_ours, errs_fred = Float64[], Float64[]
    for (a, b) in sort!(collect(keys(ours)))
        t = distance(curve(a), curve(b), ours[(a, b)])
        eo = (ours[(a, b)] - t) / t
        push!(errs_ours, eo)
        line = @sprintf("%-14s %-14s true %.12g  ours %+.2e", a, b, t, eo)
        if fred !== nothing
            ef = (fred[(a, b)] - t) / t
            push!(errs_fred, ef)
            line *= @sprintf("  Fred %+.2e", ef)
            abs(ef) > REPORT && println(line)
        end
        abs(eo) > REPORT &&
            (fred === nothing || abs(errs_fred[end]) <= REPORT) &&
            println(line)
    end
    summarize("ours", errs_ours)
    fred === nothing || summarize("Fred", errs_fred)
end

main(ARGS)
