# Exact decision procedure for the continuous Fréchet distance: the free-space
# reachability test of Alt and Godau (1995), in O(n m d) time and O(m) memory.

const FREE_EMPTY = (1.0, 0.0)

@inline free_is_empty(iv::Tuple{Float64,Float64}) = iv[1] > iv[2]

# `iv` with its lower end raised to `lo`.
@inline free_clip(iv::Tuple{Float64,Float64}, lo::Float64) = (max(iv[1], lo), iv[2])

"""
    free_interval(c, a, b, r_sq) -> (lo, hi)

The parameters `t` in [0, 1] of the points `a + t (b - a)` within distance `sqrt(r_sq)`
of `c`, as an interval, empty when `lo > hi` (`FREE_EMPTY`).

The ends are decided by the vertices where that matters: `lo == 0` exactly when
`|a - c|^2 <= r_sq`, and `hi == 1` exactly when `|b - c|^2 <= r_sq`, so whether the
free space continues through a vertex does not depend on the rounding of the roots.
The roots come from the closest point on the line (`t0`, and the squared distance to it
from a second pass), not the quadratic formula, which loses half of the significant
digits when the interval is nearly empty.
"""
@inline function free_interval(
    c::AbstractVector{T},
    a::AbstractVector{T},
    b::AbstractVector{T},
    r_sq::T,
)::Tuple{Float64,Float64} where {T}
    n = vec_length(c, a, b)
    dd = de = ee = ff = zero(T)
    # In bounds: c, a and b have n coordinates, indexed from 1 (vec_length).
    @inbounds @simd for k in 1:n
        d = b[k] - a[k]
        e = a[k] - c[k]
        f = b[k] - c[k]
        dd += d * d
        de += d * e
        ee += e * e
        ff += f * f
    end
    ends = (ee <= r_sq, ff <= r_sq)
    trivial = free_trivial(dd, ends)
    trivial === nothing || return trivial
    t0 = -de / dd
    h = zero(T)
    @inbounds @simd for k in 1:n
        g = (a[k] - c[k]) + t0 * (b[k] - a[k])
        h += g * g
    end
    return free_from_line(t0, h, dd, ends, r_sq)
end

# The free interval when it does not depend on the closest point: a segment of length 0,
# or both ends within the radius (the disk is convex); `nothing` otherwise. `ends` says
# whether the ends a and b are within the radius.
@inline function free_trivial(dd, (in_a, in_b))
    dd == 0 && return in_a ? (0.0, 1.0) : FREE_EMPTY
    (in_a && in_b) && return (0.0, 1.0)
    return nothing
end

# The free interval from the closest point `a + t0 (b - a)` on the line, the squared
# distance `h` to it, and the squared length `dd` of the segment.
@inline function free_from_line(t0::T, h::T, dd::T, (in_a, in_b), r_sq::T) where {T}
    gap = r_sq - h
    if gap < 0
        # Rounding: an endpoint inside the disk makes the interval non-empty.
        (in_a || in_b) || return FREE_EMPTY
        gap = zero(T)
    end
    half = sqrt(gap / dd)
    lo = in_a ? 0.0 : max(t0 - half, nextfloat(0.0))
    hi = in_b ? 1.0 : min(t0 + half, prevfloat(1.0))
    return (Float64(lo), Float64(hi))
end

# Where the decision reads the free space from. From the coordinates (`PointSpace`)
# every free interval costs O(d); `frechet_c_compute_exact` passes tables computed once
# (`ExactTables`), O(1) each. A source answers:
#   space_size(s) = (n, m), the numbers of vertices;
#   free_p_on_q(s, i, j, r_sq): the free interval of vertex i of P on segment j of Q;
#   free_q_on_p(s, j, i, r_sq): the free interval of vertex j of Q on segment i of P;
#   within(s, i, j, r_sq): whether |p_i - q_j|^2 <= r_sq.
struct PointSpace{V}
    p::V
    q::V
end

space_size(s::PointSpace) = (length(s.p), length(s.q))
free_p_on_q(s::PointSpace, i, j, r_sq) = free_interval(s.p[i], s.q[j], s.q[j+1], r_sq)
free_q_on_p(s::PointSpace, j, i, r_sq) = free_interval(s.q[j], s.p[i], s.p[i+1], r_sq)
within(s::PointSpace, i, j, r_sq) = fp_dist_sq(s.p[i], s.q[j]) <= r_sq

# Reachable parts of the left edges of the first column of cells (p1 against the
# segments of q): reachable from the start, and on to the next edge only through q[j+1].
function decide_first_column!(left, s, r_sq)
    reach = true
    for j in eachindex(left)
        if reach
            iv = free_p_on_q(s, 1, j, r_sq)
            left[j] = iv
            reach = iv[2] == 1.0
        else
            left[j] = FREE_EMPTY
        end
    end
    return nothing
end

# One column of cells, segment i of P against every segment of Q. `left[j]` holds the
# reachable part of the left edge of cell j and is replaced by that of its right edge;
# `bottom` is the reachable part of the bottom edge of the first cell. Returns the
# reachable part of the top edge of the last cell, and whether any edge is reachable.
function decide_column!(left, bottom, s, i::Int, r_sq)
    any_reach = !free_is_empty(bottom)
    for j in eachindex(left)
        l = left[j]
        if free_is_empty(l) && free_is_empty(bottom)
            continue  # nothing enters this cell, so nothing leaves it
        end
        top = free_q_on_p(s, j + 1, i, r_sq)
        right = free_p_on_q(s, i + 1, j, r_sq)
        # The free space of a cell is convex: entering from the left reaches the whole
        # free part of the top edge; entering only from the bottom, the part to its right.
        left[j] = free_is_empty(bottom) ? free_clip(right, l[1]) : right
        bottom = free_is_empty(l) ? free_clip(top, bottom[1]) : top
        any_reach |= !free_is_empty(left[j])
    end
    return bottom, any_reach || !free_is_empty(bottom)
end

"""
    frechet_c_decide_exact(P, Q, r) -> Bool

Exact decision procedure for the continuous Fréchet distance: is
``d_F(P, Q) \\le r``? It is the free-space reachability test of Alt and Godau (1995),
in O(n m d) time for curves of n and m vertices in ℝ^d, and O(m) memory. Unlike
`frechet_c_compute`, no tolerance is involved: the answer is that of the exact
algorithm, carried out in floating point.

Floating point can change the answer only when `r` is within rounding of the
distance (a relative 1e-15 or so, more for nearly degenerate configurations). The free
intervals are computed from the closest point on each segment, and whether the free
space continues through a vertex is decided by comparing that vertex's squared
distance with `r^2` directly.

`frechet_c_compute_exact` uses it to find the exact distance among the critical values.
Both are slow (O(n m d) per decision); `frechet_c_compute` is usually much faster.

Takes any curves the other algorithms take (static or `Vector` points; matrix curves
are converted). Curves of one vertex are allowed. `r` must be finite and non-negative.

# Example
```julia
m = frechet_c_compute(P, Q)
frechet_c_decide_exact(P, Q, m.leash)                 # true
frechet_c_decide_exact(P, Q, m.leash * (1 - 1e-6))    # false, unless m.leash is not optimal
```
"""
function frechet_c_decide_exact(
    P::GPntSeq{N,T,PT},
    Q::GPntSeq{N,T,PT},
    r::Real,
)::Bool where {N,T,PT<:AbstractVector{T}}
    check_curves(P, Q, "frechet_c_decide_exact"; min_vertices = 1)
    check_nonnegative(r, "frechet_c_decide_exact", "the distance r")
    return decide_free_space(P.pnts, Q.pnts, T(r)^2)
end

# The decision of frechet_c_decide_exact on vertex vectors p and q, for the squared
# radius r_sq, without the input checks.
decide_free_space(p, q, r_sq::Float64)::Bool = decide_free_space(PointSpace(p, q), r_sq)

# The decision on a free-space source (see PointSpace): O(n m) free intervals.
function decide_free_space(s, r_sq::Float64)::Bool
    n, m = space_size(s)
    (within(s, 1, 1, r_sq) && within(s, n, m, r_sq)) || return false
    n == 1 && return all(j -> within(s, 1, j, r_sq), 1:m)
    m == 1 && return all(i -> within(s, i, 1, r_sq), 1:n)

    left = Vector{Tuple{Float64,Float64}}(undef, m - 1)
    decide_first_column!(left, s, r_sq)
    open_row = true  # the bottom edges of the first row are reachable so far
    top = FREE_EMPTY
    for i in 1:(n-1)
        bottom = FREE_EMPTY
        if open_row
            bottom = free_q_on_p(s, 1, i, r_sq)
            open_row = bottom[2] == 1.0
        end
        top, any_reach = decide_column!(left, bottom, s, i, r_sq)
        any_reach || return false
    end
    # The end (n, m) is reachable through the last cell's right edge or its top edge.
    return (!free_is_empty(left[m-1]) && left[m-1][2] == 1.0) ||
           (!free_is_empty(top) && top[2] == 1.0)
end
