# Exact decision procedure for the continuous Fréchet distance: the free-space
# reachability test of Alt and Godau (1995), in O(n m d) time and O(n + m) memory,
# with exact shortcuts that settle most queries much faster.

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
Base.@propagate_inbounds free_p_on_q(s::PointSpace, i, j, r_sq) =
    free_interval(s.p[i], s.q[j], s.q[j+1], r_sq)
Base.@propagate_inbounds free_q_on_p(s::PointSpace, j, i, r_sq) =
    free_interval(s.q[j], s.p[i], s.p[i+1], r_sq)
within(s::PointSpace, i, j, r_sq) = fp_dist_sq(s.p[i], s.q[j]) <= r_sq

# Reachable parts of the left edges of the first column of cells (p1 against the
# segments of q): reachable from the start, and on to the next edge only through q[j+1].
function decide_first_column!(left, s, r_sq)
    reach = true
    for j in eachindex(left)
        if reach
            # In bounds: j is an edge of `left`, so q[j + 1] exists.
            iv = @inbounds free_p_on_q(s, 1, j, r_sq)
            left[j] = iv
            reach = iv[2] == 1.0
        else
            left[j] = FREE_EMPTY
        end
    end
    return nothing
end

# A free-space source and the squared radius of one decision.
struct FreeSweep{S}
    s::S
    r_sq::Float64
end

# The reachable parts of the right and top edges of cell (i, j), entered through the
# reachable parts l of its left edge and b of its bottom edge. The free space of a cell
# is convex: entering from the left reaches the whole free part of the top edge;
# entering only from the bottom, the part to its right (and the same for the right
# edge, the other way around).
Base.@propagate_inbounds function cell_exits(sw::FreeSweep, i::Int, j::Int, l, b)
    top = free_q_on_p(sw.s, j + 1, i, sw.r_sq)
    right = free_p_on_q(sw.s, i + 1, j, sw.r_sq)
    return (free_is_empty(b) ? free_clip(right, l[1]) : right),
    (free_is_empty(l) ? free_clip(top, b[1]) : top)
end

# One column of cells, segment i of P against the segments of Q. `left[j]` holds the
# reachable part of the left edge of cell j for j in `band` (the edges outside it are
# unreachable), and is replaced by that of its right edge; `bottom` is the reachable part
# of the bottom edge of the first cell. Only the cells that something enters are visited:
# from the first reachable left edge (or the first cell, when `bottom` is reachable) up
# to the last cell reached from the left or from below. Returns the reachable part of the
# top edge of the last cell, whether any edge is reachable, and the band of the next
# column (empty when no right edge is reachable).
function decide_column!(left, bottom, sw::FreeSweep, i::Int, band::UnitRange{Int})
    any_reach = !free_is_empty(bottom)
    lo = hi = 0  # the reachable right edges, lo:hi
    j = free_is_empty(bottom) ? first(band) : 1
    while j <= length(left)
        l = j in band ? left[j] : FREE_EMPTY
        if free_is_empty(l) && free_is_empty(bottom)
            j > last(band) && break  # nothing enters this cell or any cell above it
            left[j] = FREE_EMPTY  # nothing enters this cell, so nothing leaves it
        else
            # In bounds: j <= m - 1 and i <= n - 1, so q[j + 1] and p[i + 1] exist.
            left[j], bottom = @inbounds cell_exits(sw, i, j, l, bottom)
            free_is_empty(left[j]) || ((lo, hi) = (lo == 0 ? j : lo, j))
        end
        j += 1
    end
    # After a break the top edge of the last cell is unreachable (bottom is empty).
    any_reach |= lo != 0
    return bottom, any_reach || !free_is_empty(bottom), lo == 0 ? (1:0) : (lo:hi)
end

# Whether the end (n, m) is reachable, through the right edge of the last cell (the
# last of `left`, if in `band`) or its top edge `top`.
function end_reached(left, band::UnitRange{Int}, top)::Bool
    k = length(left)
    return (k in band && !free_is_empty(left[k]) && left[k][2] == 1.0) ||
           (!free_is_empty(top) && top[2] == 1.0)
end

# The hierarchies of curve pieces of the free-space search: their bounding boxes (BoxTree
# in frechet_decide_blocks.jl, LazyBoxTree in frechet_decide_lazy_tree.jl), or sausages
# (SausageTree in frechet_decide_sausage.jl).
abstract type PieceTree end
abstract type AbstractBoxTree <: PieceTree end

"""
    frechet_c_decide_exact(P, Q, r; tree = nothing) -> Bool

Exact decision procedure for the continuous Fréchet distance: is
``d_F(P, Q) \\le r``? It is the free-space reachability test of Alt and Godau (1995),
in O(n m d) time for curves of n and m vertices in ℝ^d, and O(n + m) memory. Unlike
`frechet_c_compute`, no tolerance is involved: the answer is that of the exact
algorithm, carried out in floating point.

Floating point can change the answer only when `r` is within rounding of the
distance (a relative 1e-15 or so, more for nearly degenerate configurations). The free
intervals are computed from the closest point on each segment, and whether the free
space continues through a vertex is decided by comparing that vertex's squared
distance with `r^2` directly.

Most queries are settled sooner, by tests that are exact as well:

- the endpoint distances, and an upper bound (the endpoint distance plus the widths of
  the curves);
- greedy monotone walks through pairs of vertices, one step at a time or by adaptive
  steps over pieces of the curves that lie in balls within `r` of each other, and, in
  dimension 6 and up, over pieces close (in the Fréchet distance) to segments whose
  ends are within `r` (a "yes");
- a vertex of one curve with no point of the other within `r`, found by a descent of a
  tree of the pieces of the curve (a "no");
- decisions on simplifications of the curves, which settle the query when the
  simplified curves are far enough from `r`;
- in the free-space search, a recursive decomposition of the diagram into boxes:
  a box whose two curve pieces are all within `r` of each other is settled through its
  boundary, and a box that nothing enters is skipped.

The pieces of the curves are described by `tree`: `LazyBoxTree` (their bounding boxes,
computed as the search needs them), `BoxTree` (the same, all computed first), or
`SausageTree` (the segment between the ends of a piece, thickened to contain it: a
description that needs only distances between points, so it does not grow with the
dimension). The default, `nothing`, takes `LazyBoxTree` below dimension
`SAUSAGE_TREE_DIM` (64) and `SausageTree` from it on (see `default_tree`).

The answer is that of the plain sweep, up to rounding at the distance itself. The
worst case, a query close to the distance between curves that wander within `r` of
each other, is still O(n m d).

`frechet_c_compute_exact` uses it to find the exact distance among the critical values.

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
    r::Real;
    tree::Union{Nothing,Type{<:PieceTree}} = nothing,
)::Bool where {N,T,PT<:AbstractVector{T}}
    check_curves(P, Q, "frechet_c_decide_exact"; min_vertices = 1)
    check_nonnegative(r, "frechet_c_decide_exact", "the distance r")
    # The ends must be matched (the sweep's own test, in O(1)).
    r_sq = Float64(r)^2
    (fp_dist_sq(P[1], Q[1]) <= r_sq && fp_dist_sq(P[end], Q[end]) <= r_sq) || return false
    # d_F(P, Q) is at most the endpoint distance plus the widths of the curves (O(n + m)).
    # The margin covers the rounding of that bound, so a "yes" here is never wrong.
    if (length(P) >= 2) &&
       (length(Q) >= 2) &&
       frechet_dist_upper_bound(P, Q) * (1 + 1e-12) <= r
        return true
    end
    kind = tree === nothing ? default_tree(P.pnts) : tree
    return decide_simplified(P, Q, Float64(r), kind)
end

"""
    DecideCurve(P, Tree = default_tree(P.pnts))

A curve prepared for many calls of `frechet_c_decide_exact`: it is checked once, and
its width (`frechet_width_approx`) and bounding box are stored, so that the O(1) filters
of the decision do not touch its vertices. Its prefix lengths and the bounding boxes of a
hierarchy of its pieces (for the free-space search) are computed when a decision first
needs them, and kept: most decisions need little or none of the hierarchy. Preparing
takes O(n d) time, and the curve O(n d) memory.

`Tree` is the description of the pieces (see `frechet_c_decide_exact`): by default
`LazyBoxTree` below dimension `SAUSAGE_TREE_DIM` and `SausageTree` from it on. Both are
updated by the decisions that use them, so such a DecideCurve must not be used by several
threads at once. `Tree = BoxTree` computes everything when the curve is prepared; such a
DecideCurve is never written to, and may be shared by threads.

`frechet_c_decide_exact(A::DecideCurve, B::DecideCurve, r)` gives the same answer as
`frechet_c_decide_exact(A.P, B.P, r)`.
"""
struct DecideCurve{C<:GPntSeq,TR}
    P::C
    width::Float64
    lo::Vector{Float64}
    hi::Vector{Float64}
    tree::TR  # the tree of the segments of P, for the free-space search
end

DecideCurve(P::GPntSeq) = DecideCurve(P, default_tree(P.pnts))

function DecideCurve(
    P::GPntSeq{N,T,PT},
    ::Type{TR},
) where {N,T,PT<:AbstractVector{T},TR<:PieceTree}
    check_curve(P, "DecideCurve", "the curve"; min_vertices = 1)
    d = length(P[1])
    lo = fill(Inf, d)
    hi = fill(-Inf, d)
    for x in P, c in 1:d
        lo[c] = min(lo[c], x[c])
        hi[c] = max(hi[c], x[c])
    end
    w = length(P) >= 2 ? Float64(frechet_width_approx(P)) : 0.0
    return DecideCurve(P, w, lo, hi, TR(P.pnts))
end

# The largest squared distance between a point of A's box and a point of B's.
function box_far_sq(A::DecideCurve, B::DecideCurve)::Float64
    far = 0.0
    for c in eachindex(A.lo)
        span = max(A.hi[c] - B.lo[c], B.hi[c] - A.lo[c])
        far += span * span
    end
    return far
end

function frechet_c_decide_exact(
    A::DecideCurve{C},
    B::DecideCurve{C},
    r::Real,
)::Bool where {C}
    check_nonnegative(r, "frechet_c_decide_exact", "the distance r")
    length(A.lo) == length(B.lo) || throw(
        ArgumentError(
            "frechet_c_decide_exact: the two curves must have the same dimension",
        ),
    )
    P, Q = A.P, B.P
    r_sq = Float64(r)^2
    e1, e2 = fp_dist_sq(P[1], Q[1]), fp_dist_sq(P[end], Q[end])
    (e1 <= r_sq && e2 <= r_sq) || return false
    # Every point of one curve within r of every point of the other: any matching works.
    box_far_sq(A, B) * (1 + 1e-9) <= r_sq && return true
    if (length(P) >= 2) &&
       (length(Q) >= 2) &&
       (A.width + B.width + sqrt(max(e1, e2))) * (1 + 1e-12) <= r
        return true
    end
    return decide_simplified(P, Q, Float64(r), (A.tree, B.tree))
end

# The decision through simplifications, coarse first. simplify_vertices(p, w^2) keeps a
# vertex when it is farther than w from the last kept one, so d_F(P, P') <= w: the skipped
# vertices and the segments between them lie within w of the last kept vertex, which they
# are matched to. Hence d_F(P', Q') <= r - 2w implies d_F(P, Q) < r, and
# d_F(P', Q') > r + 2w implies d_F(P, Q) > r. The margin on 2w covers the rounding of the
# distances to the kept vertices. Rounds stop once a simplification keeps more than a
# quarter of a curve; then the curves themselves are decided.
# `trees` are the trees of P and Q (see PieceTree), or the kind of tree to build;
# the simplifications get trees of the same kind.
function decide_simplified(P, Q, r::Float64, trees = LazyBoxTree)::Bool
    p, q = P.pnts, Q.pnts
    kind = trees isa Tuple ? tree_kind(trees[1]) : trees
    walked = false
    if length(p) > 2 && length(q) > 2
        tp, tq = tree_pair(p, q, trees)
        # The greedy walks prove most "yes" answers; then a vertex with no point of the
        # other curve within r proves many "no" answers; then, in higher dimensions, the
        # walk by sausages, whose steps cost scans of the pieces, proves more "yes" answers.
        greedy_vertex_walk(p, q, r^2) && return true
        greedy_ball_walk(p, q, tree_lens(tp), tree_lens(tq), r) && return true
        (far_vertex(p, q, tq, r^2) || far_vertex(q, p, tp, r^2)) && return false
        box_dim(p) >= SAUSAGE_WALK_DIM &&
            greedy_sausage_walk(p, q, r, tp, tq) &&
            return true
        trees = (tp, tq)
        walked = true
    end
    w = r / 8.0
    while w > 0.0
        ps = simplify_vertices(p, w * w, length(p) ÷ 8)
        ps === nothing && break
        qs = simplify_vertices(q, w * w, length(q) ÷ 8)
        qs === nothing && break
        e = 2.0 * w * (1.0 + 1.0e-9)
        decide_free_space(ps, qs, (r + e)^2, kind) || return false
        (r > e) && decide_free_space(ps, qs, (r - e)^2, kind) && return true
        w /= 8.0
    end
    return decide_free_space(p, q, r^2, trees, walked)
end

# The vertices p[1], the ones whose squared distance to the last kept vertex exceeds w_sq,
# and p[end]; `nothing` as soon as more than `limit` (at least 2) would be kept.
function simplify_vertices(p::AbstractVector, w_sq::Float64, limit::Int)
    limit >= 2 || return nothing
    out = [p[1]]
    cur = p[1]
    for i in 2:(length(p)-1)
        fp_dist_sq(p[i], cur) > w_sq || continue
        cur = p[i]
        push!(out, cur)
        length(out) >= limit && return nothing  # with p[end], more than limit
    end
    push!(out, p[end])
    return out
end

# The decision of frechet_c_decide_exact on vertex vectors p and q, for the squared
# radius r_sq, without the input checks. The greedy walk answers many "yes" instances in
# O(n + m); the free space is searched only when it fails.
# `trees`: as for decide_simplified; `walked`: the greedy walks have been tried already.
function decide_free_space(
    p,
    q,
    r_sq::Float64,
    trees = LazyBoxTree,
    walked::Bool = false,
)::Bool
    !walked && greedy_vertex_walk(p, q, r_sq) && return true
    s = PointSpace(p, q)
    n, m = space_size(s)
    if n - 1 > DECIDE_BLOCK && m - 1 > DECIDE_BLOCK
        (within(s, 1, 1, r_sq) && within(s, n, m, r_sq)) || return false
        tp, tq = tree_pair(p, q, trees)
        !walked &&
            greedy_ball_walk(p, q, tree_lens(tp), tree_lens(tq), sqrt(r_sq)) &&
            return true
        return decide_free_space_blocks(p, q, r_sq, DECIDE_BLOCK, (tp, tq))
    end
    return decide_free_space(s, r_sq)
end

# The decision on a free-space source (see PointSpace): O(n m) free intervals.
function decide_free_space(s, r_sq::Float64)::Bool
    n, m = space_size(s)
    (within(s, 1, 1, r_sq) && within(s, n, m, r_sq)) || return false
    n == 1 && return all(j -> within(s, 1, j, r_sq), 1:m)
    m == 1 && return all(i -> within(s, i, 1, r_sq), 1:n)

    left = Vector{Tuple{Float64,Float64}}(undef, m - 1)
    decide_first_column!(left, s, r_sq)
    sw = FreeSweep(s, r_sq)
    band = 1:(m-1)  # the left edges that may be reachable
    open_row = true  # the bottom edges of the first row are reachable so far
    top = FREE_EMPTY
    for i in 1:(n-1)
        bottom = FREE_EMPTY
        if open_row
            # In bounds: i <= n - 1, so p[i + 1] exists.
            bottom = @inbounds free_q_on_p(s, 1, i, r_sq)
            open_row = bottom[2] == 1.0
        end
        top, any_reach, band = decide_column!(left, bottom, sw, i, band)
        any_reach || return false
    end
    return end_reached(left, band, top)
end

# Whether a monotone path through pairs of vertices stays within the radius, found
# greedily (to the closest of the next pairs). Such a path gives d_F(P, Q) <= r exactly:
# along a step the matched points move linearly, so their distance is convex and at
# most that of the two pairs at its ends. O(n + m) distances; false says nothing.
function greedy_vertex_walk(p, q, r_sq::Float64)::Bool
    n, m = length(p), length(q)
    i = j = 1
    fp_dist_sq(p[1], q[1]) <= r_sq || return false
    while (i < n) || (j < m)
        d_both = (i < n && j < m) ? fp_dist_sq(p[i+1], q[j+1]) : Inf
        d_p = i < n ? fp_dist_sq(p[i+1], q[j]) : Inf
        d_q = j < m ? fp_dist_sq(p[i], q[j+1]) : Inf
        d = min(d_both, d_p, d_q)
        d <= r_sq || return false
        if d == d_both
            i += 1
            j += 1
        elseif d == d_p
            i += 1
        else
            j += 1
        end
    end
    return true
end
