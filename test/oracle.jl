# Independent reference implementations used to test the Frechet algorithms.
# They share no code with FrechetDist: points are plain coordinate vectors.

using LinearAlgebra: dot

oracle_dist(p, q) = sqrt(sum(abs2, p .- q))

"""
    oracle_discrete(P, Q; diagonal = true)

Discrete Frechet distance of two vectors of points, by the classic dynamic program.
With `diagonal = false` only single steps are allowed (a "staircase" coupling), which
is the variant FrechetDist's discrete routines compute.
"""
function oracle_discrete(P, Q; diagonal = true)
    n, m = length(P), length(Q)
    c = fill(Inf, n, m)
    for i = 1:n, j = 1:m
        best = 0.0
        if i > 1 || j > 1
            best = min(
                i > 1 ? c[i-1, j] : Inf,
                j > 1 ? c[i, j-1] : Inf,
                (diagonal && i > 1 && j > 1) ? c[i-1, j-1] : Inf,
            )
        end
        c[i, j] = max(oracle_dist(P[i], Q[j]), best)
    end
    return c[n, m]
end

"""
    oracle_dtw(P, Q; diagonal = true)

Discrete dynamic time warping price: the cheapest monotone vertex matching, where a
matching costs the sum of the distances of its matched pairs. `diagonal = false`
restricts the matching to single steps, as FrechetDist's `DTW_d_compute` does.
"""
function oracle_dtw(P, Q; diagonal = true)
    n, m = length(P), length(Q)
    c = fill(Inf, n, m)
    for i = 1:n, j = 1:m
        best = 0.0
        if i > 1 || j > 1
            best = min(
                i > 1 ? c[i-1, j] : Inf,
                j > 1 ? c[i, j-1] : Inf,
                (diagonal && i > 1 && j > 1) ? c[i-1, j-1] : Inf,
            )
        end
        c[i, j] = oracle_dist(P[i], Q[j]) + best
    end
    return c[n, m]
end

"""
Interval of t in [0, 1] with |a + t (b - a) - c| <= eps, or `nothing`.

Computed from the closest point on the line instead of the quadratic formula, which
loses half of the significant digits when the interval is nearly empty.
"""
function free_interval(a, b, c, eps)
    d = b .- a
    f = a .- c
    A = dot(d, d)
    if A == 0
        return dot(f, f) <= eps^2 ? (0.0, 1.0) : nothing
    end
    t0 = -dot(f, d) / A
    r = f .+ t0 .* d
    gap = eps^2 - dot(r, r)
    gap < 0 && return nothing
    half = sqrt(gap / A)
    lo = max(t0 - half, 0.0)
    hi = min(t0 + half, 1.0)
    return lo <= hi ? (lo, hi) : nothing
end

# Reachable part of the next edge when crossing a cell; `entry` is the reachable
# part of the edge on the near side, `other` the reachable part of the other entry
# edge, `free` the free interval of the exit edge.
function exit_interval(entry, other, free)
    free === nothing && return nothing
    entry !== nothing && return free
    other === nothing && return nothing
    lo = max(free[1], other[1])
    return lo <= free[2] ? (lo, free[2]) : nothing
end

"Is the continuous Frechet distance at most eps? (Alt-Godau free-space decision.)"
function oracle_decide(P, Q, eps)
    n, m = length(P), length(Q)
    oracle_dist(P[1], Q[1]) <= eps && oracle_dist(P[n], Q[m]) <= eps || return false
    (n == 1 || m == 1) &&
        return all(p -> oracle_dist(p, n == 1 ? P[1] : Q[1]) <= eps, n == 1 ? Q : P)

    # LF[i][j]: free part of segment Q[j]Q[j+1] seen from P[i]; BF[i][j]: free part of
    # segment P[i]P[i+1] seen from Q[j].
    LF = [free_interval(Q[j], Q[j+1], P[i], eps) for i = 1:n, j = 1:(m-1)]
    BF = [free_interval(P[i], P[i+1], Q[j], eps) for i = 1:(n-1), j = 1:m]
    LR = Matrix{Any}(nothing, n, m - 1)   # reachable parts of the left edges
    BR = Matrix{Any}(nothing, n - 1, m)   # reachable parts of the bottom edges
    LR[1, 1] = LF[1, 1]
    BR[1, 1] = BF[1, 1]
    for j = 2:(m-1)
        prev = LR[1, j-1]
        LR[1, j] =
            (
                prev !== nothing &&
                prev[2] == 1.0 &&
                LF[1, j] !== nothing &&
                LF[1, j][1] == 0.0
            ) ? LF[1, j] : nothing
    end
    for i = 2:(n-1)
        prev = BR[i-1, 1]
        BR[i, 1] =
            (
                prev !== nothing &&
                prev[2] == 1.0 &&
                BF[i, 1] !== nothing &&
                BF[i, 1][1] == 0.0
            ) ? BF[i, 1] : nothing
    end
    for i = 1:(n-1), j = 1:(m-1)
        left, bottom = LR[i, j], BR[i, j]
        LR[i+1, j] = exit_interval(bottom, left, LF[i+1, j])
        BR[i, j+1] = exit_interval(left, bottom, BF[i, j+1])
    end
    final = LR[n, m-1]
    return final !== nothing && final[2] == 1.0
end

"""
    oracle_continuous(P, Q; rtol = 1e-12)

Continuous Frechet distance of two vectors of points: bisection on the free-space
decision procedure, between the end-to-end distances and the discrete distance.
"""
function oracle_continuous(P, Q; rtol = 1e-12)
    lo = max(oracle_dist(P[1], Q[1]), oracle_dist(P[end], Q[end]))
    hi = oracle_discrete(P, Q) * (1 + 1e-12)
    oracle_decide(P, Q, lo) && return lo
    while hi - lo > rtol * hi
        mid = (lo + hi) / 2
        oracle_decide(P, Q, mid) ? (hi = mid) : (lo = mid)
    end
    return hi
end
