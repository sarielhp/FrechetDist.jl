# The free-space decision of frechet_decide.jl, by a recursive decomposition of the
# diagram into boxes. A box is a range of cells (segments of P) by a range of rows
# (segments of Q):
#
# - a box that nothing enters is skipped (its exits are as empty as its entries);
# - a box whose two curve pieces are within r of each other (every point of one within r
#   of every point of the other, by their bounding boxes, or their sausages in
#   frechet_decide_sausage.jl) is free, and what is reachable on its top and right edges
#   follows from where it is entered (the free space of the box is convex);
# - otherwise it is split in four, or in two across its longer side, down to boxes of at
#   most DECIDE_BLOCK x DECIDE_BLOCK cells, whose cells are visited one by one, as in
#   decide_column!.
#
# (A box that is entered has a free point on its boundary, so it is never empty: there is
# no test for boxes whose pieces are apart.) Free boxes take time proportional to their
# perimeter, so large uniform parts of the free space cost little. The ranges are those
# of a segment tree over the segments of each curve, whose nodes describe their pieces.
#
# The test has a relative margin of 1e-9, so a box is free only when the cell by cell
# computation would give every free interval (0, 1) as well: the answers are those of
# decide_free_space(::PointSpace, r_sq).
#
# Throughout, `left[j]` holds the reachable part of the left edge of the first
# unprocessed cell of row j, and `bot[i]` that of the bottom edge of the first
# unprocessed cell of column i. Processing a box replaces its left and bottom edges by
# its right and top edges, so the two halves of a split box are processed in order (left
# then right, or bottom then top) on the same two arrays.

"""Largest side of the boxes of the free-space decision whose cells are visited one by one;
also the number of segments of the bottom nodes of the trees (BoxTree, LazyBoxTree)."""
const DECIDE_BLOCK = 8

# A segment tree over the segments 1:ncell of a curve: its N bottom nodes (N a power of
# two) are the runs of DECIDE_BLOCK consecutive segments, node 1 is all of them, the
# children of node v are 2v and 2v + 1, and its parent v >> 1 (so the tree needs no
# pointers). Column v of
# box is the bounding box of the piece of the curve of the range of node v, in one block:
# its d smallest coordinates, then its d largest (Inf and -Inf past ncell).
# lens[k] is the length of the curve up to vertex k. All of it is computed when the tree is
# built; LazyBoxTree (frechet_decide_lazy_tree.jl) computes the same boxes on first use.
# Readers call node_box!(t, v) before reading column v, and tree_lens(t) for lens.
struct BoxTree <: AbstractBoxTree
    box::Matrix{Float64}  # 2d x 2N
    N::Int
    ncell::Int
    lens::Vector{Float64}
end

# The dimension of the points of p: a compile-time constant for static points, so that
# the loops over coordinates unroll.
box_dim(p::AbstractVector{<:StaticVector{D}}) where {D} = D
box_dim(p) = length(p[1])

function BoxTree(p)
    d = box_dim(p)
    ncell = length(p) - 1
    N = nextpow(2, max(cld(ncell, DECIDE_BLOCK), 1))
    box = Matrix{Float64}(undef, 2d, 2N)
    t = BoxTree(box, N, ncell, zeros(ncell + 1))
    for w in N:(2N-1)
        scan_box!(box, p, node_range(t, w), w)
    end
    for v in (N-1):-1:1
        union_box!(box, v, d)
    end
    lens = t.lens
    for c in 1:ncell
        lens[c+1] = lens[c] + fp_dist(p[c], p[c+1])
    end
    return t
end

# Column w of box: the bounding box of the vertices of the segments rng of p, by one pass
# over them (Inf and -Inf when rng is empty). For static points the box is kept in two
# points.
function scan_box!(
    box::Matrix{Float64},
    p::AbstractVector{<:StaticVector{D}},
    rng::UnitRange{Int},
    w::Int,
) where {D}
    if isempty(rng)
        lo, hi = fill(Inf, p[1]), fill(-Inf, p[1])
    else
        # In bounds: the vertices of the segments rng are within p.
        lo = hi = @inbounds p[first(rng)]
        @inbounds for c in (first(rng)+1):(last(rng)+1)
            x = p[c]
            lo = box_min.(lo, x)
            hi = box_max.(hi, x)
        end
    end
    # In bounds: w is a column of box, which has 2D rows.
    @inbounds for k in 1:D
        box[k, w], box[D+k, w] = lo[k], hi[k]
    end
    return nothing
end

function scan_box!(box::Matrix{Float64}, p, rng::UnitRange{Int}, w::Int)
    d = box_dim(p)
    verts = isempty(rng) ? (1:0) : (first(rng):(last(rng)+1))
    # In bounds: the vertices of the segments rng are within p, and w is a column of box.
    @inbounds for k in 1:d
        lo, hi = Inf, -Inf
        for c in verts
            x = p[c][k]
            lo = box_min(lo, x)
            hi = box_max(hi, x)
        end
        box[k, w], box[d+k, w] = lo, hi
    end
    return nothing
end

# The box of node v from those of its children (d: the dimension).
@inline function union_box!(box::Matrix{Float64}, v::Int, d::Int)
    # In bounds: v has children, so 2v + 1 is a column of box.
    @inbounds for k in 1:d
        box[k, v] = box_min(box[k, 2v], box[k, 2v+1])
        box[d+k, v] = box_max(box[d+k, 2v], box[d+k, 2v+1])
    end
    return nothing
end

node_box!(::BoxTree, ::Int) = nothing

# min and max of finite coordinates (the curves are checked), without the NaN and signed
# zero handling of Base.min and Base.max.
box_min(a::T, b::T) where {T<:Real} = ifelse(b < a, b, a)
box_max(a::T, b::T) where {T<:Real} = ifelse(b > a, b, a)
tree_lens(t::BoxTree) = t.lens
tree_kind(::BoxTree) = BoxTree

# The trees of p and q: given, or built.
tree_pair(p, q, trees::Tuple) = trees
tree_pair(p, q, ::Type{T}) where {T<:PieceTree} = (T(p), T(q))

# The number of cells a node of depth-index v spans (before clipping at the end).
node_size(t::PieceTree, v::Int) = (t.N >> (63 - leading_zeros(v))) * DECIDE_BLOCK

# Whether v is a bottom node (it has no children).
is_bottom(t::PieceTree, v::Int) = v >= t.N

# The cells of node v (possibly empty, past the end of the curve).
function node_range(t::PieceTree, v::Int)::UnitRange{Int}
    h = 63 - leading_zeros(v)  # the depth of v
    size = (t.N >> h) * DECIDE_BLOCK
    first = (v - (1 << h)) * size + 1
    return first:min(first+size-1, t.ncell)
end

# Whether every point of the piece of node u is within the radius of every point of the
# piece of node v: by the largest squared distance between their bounding boxes.
function pieces_within(
    tp::AbstractBoxTree,
    tq::AbstractBoxTree,
    u::Int,
    v::Int,
    r_sq::Float64,
)::Bool
    node_box!(tp, u)
    node_box!(tq, v)
    far = 0.0
    A, B = tp.box, tq.box
    d = size(A, 1) >> 1
    # In bounds: u and v are nodes of their trees, whose matrices have the same rows.
    @inbounds for c in 1:d
        span = max(A[d+c, u] - B[c, v], B[d+c, v] - A[c, u])
        far += span * span
    end
    return far * (1.0 + 1.0e-9) <= r_sq
end

# One exit side of a full box, e[r]: its right edges (in `left`) or its top edges (in
# `bot`), which hold the entries through the opposite side on input. Entered through
# the adjacent side (`across`), all of it is reachable; otherwise it is reachable from
# its first entry on (entered from the left at row j, the right edges from j up).
function full_side!(e, r::UnitRange{Int}, across::Bool)
    if across
        fill!(view(e, r), (0.0, 1.0))
        return nothing
    end
    k = findfirst(t -> !free_is_empty(e[t]), r)
    k === nothing && return nothing
    t0 = r[k]
    lo = e[t0][1]
    for t in r
        e[t] = t < t0 ? FREE_EMPTY : (t == t0 ? (lo, 1.0) : (0.0, 1.0))
    end
    return nothing
end

# A full box: the reachable parts of its right edges (left[rows]) and top edges
# (bot[cols]) from where it is entered.
function full_box!(left, bot, rows::UnitRange{Int}, cols::UnitRange{Int})
    from_left = any(j -> !free_is_empty(left[j]), rows)
    from_below = any(i -> !free_is_empty(bot[i]), cols)
    full_side!(left, rows, from_below)
    full_side!(bot, cols, from_left)
    return nothing
end

# The cells (i, rows) one by one, as in decide_column!: `left` holds their left edges
# and receives their right edges; returns the top edge of the last one.
function column_cells!(left, bottom, sw::FreeSweep, i::Int, rows::UnitRange{Int})
    # In bounds: rows lies in 1:(m - 1) (the edges of `left`) and i in 1:(n - 1), so the
    # vertices i, i + 1, j and j + 1 exist.
    @inbounds for j in rows
        l = left[j]
        if free_is_empty(l) && free_is_empty(bottom)
            continue  # nothing enters this cell, so nothing leaves it
        end
        left[j], bottom = cell_exits(sw, i, j, l, bottom)
    end
    return bottom
end

# The reachable bottom edges of the first row: reachable from the start, and on to the
# next edge only through its end.
function first_row_bottoms(s, n::Int, r_sq)
    row = fill(FREE_EMPTY, n - 1)
    for i in 1:(n-1)
        # In bounds: i <= n - 1, so p[i + 1] exists.
        row[i] = @inbounds free_q_on_p(s, 1, i, r_sq)
        row[i][2] == 1.0 || break
    end
    return row
end

# Whether every vertex of p[is] is within the radius of every vertex of q[js], with a
# margin; then every free interval of the cells between them is (0, 1), the distance
# being convex along the segments.
function vertices_within(sw::FreeSweep, is::UnitRange{Int}, js::UnitRange{Int})::Bool
    p, q = sw.s.p, sw.s.q
    lim = sw.r_sq * (1.0 - 1.0e-9)
    for i in is, j in js
        fp_dist_sq(p[i], q[j]) <= lim || return false
    end
    return true
end

# The recursive decomposition of one decision.
struct BoxSweep{S,TP,TQ}
    sweep::FreeSweep{S}
    tp::TP
    tq::TQ
    leaf::Int
end

# Box (u, v) (node u of P's tree, node v of Q's), entered through left[rows] and
# bot[cols]: replaces them by its right and top edges.
function box_step!(left, bot, bs::BoxSweep, u::Int, v::Int)
    cols, rows = node_range(bs.tp, u), node_range(bs.tq, v)
    (isempty(cols) || isempty(rows)) && return nothing
    entered =
        any(i -> !free_is_empty(bot[i]), cols) || any(j -> !free_is_empty(left[j]), rows)
    entered || return nothing  # its exits are as empty as its entries
    if pieces_within(bs.tp, bs.tq, u, v, bs.sweep.r_sq)
        full_box!(left, bot, rows, cols)
    elseif (length(cols) <= bs.leaf || is_bottom(bs.tp, u)) &&
           (length(rows) <= bs.leaf || is_bottom(bs.tq, v))
        if vertices_within(bs.sweep, first(cols):(last(cols)+1), first(rows):(last(rows)+1))
            full_box!(left, bot, rows, cols)
            return nothing
        end
        for i in cols
            bot[i] = column_cells!(left, bot[i], bs.sweep, i, rows)
        end
    else
        split_box!(left, bot, bs, u, v)
    end
    return nothing
end

# Box (u, v) by its parts: in four when its sides are of the same size (bottom-left,
# bottom-right, top-left, top-right: each part after the ones left of it and below it),
# and otherwise in two across its longer side (left then right, or bottom then top).
function split_box!(left, bot, bs::BoxSweep, u::Int, v::Int)
    su, sv = node_size(bs.tp, u), node_size(bs.tq, v)
    if su == sv
        box_step!(left, bot, bs, 2u, 2v)
        box_step!(left, bot, bs, 2u + 1, 2v)
        box_step!(left, bot, bs, 2u, 2v + 1)
        box_step!(left, bot, bs, 2u + 1, 2v + 1)
    elseif su > sv
        box_step!(left, bot, bs, 2u, v)
        box_step!(left, bot, bs, 2u + 1, v)
    else
        box_step!(left, bot, bs, u, 2v)
        box_step!(left, bot, bs, u, 2v + 1)
    end
    return nothing
end

# decide_free_space(PointSpace(p, q), r_sq) by boxes; n, m > 1 and the ends checked.
# Boxes of at most B x B cells are visited cell by cell. `trees` are the trees of p and
# q, or the kind of tree to build (see tree_pair).
function decide_free_space_blocks(p, q, r_sq::Float64, B::Int, trees = LazyBoxTree)::Bool
    s = PointSpace(p, q)
    n, m = length(p), length(q)
    tp, tq = tree_pair(p, q, trees)
    bs = BoxSweep(FreeSweep(s, r_sq), tp, tq, B)
    left = Vector{Tuple{Float64,Float64}}(undef, m - 1)
    decide_first_column!(left, s, r_sq)
    bot = first_row_bottoms(s, n, r_sq)
    box_step!(left, bot, bs, 1, 1)
    return end_reached(left, 1:(m-1), bot[n-1])
end

# A ball containing the piece p[a:b]: around its middle vertex, of radius the longer of
# the two arc lengths from it (lp: the prefix lengths of p).
function piece_ball(p, lp, a::Int, b::Int)
    m = (a + b + 1) >> 1
    return p[m], max(lp[m] - lp[a], lp[b] - lp[m])
end

# Whether every point of one ball is within r of every point of the other.
function balls_within((cp, rp), (cq, rq), r::Float64)::Bool
    s = r - (rp + rq)
    return s >= 0.0 && fp_dist_sq(cp, cq) <= s * s
end

# A monotone walk from (1, 1) to (n, m) through pairs of vertices, by adaptive steps:
# from (i, j) to (i', j') when the pieces p[i:i'] and q[j:j'] are within r of each other
# (by their balls, piece_ball), the steps advancing both curves at proportional rates and growing on
# success, shrinking on failure; at the smallest step, to the closest of the next
# vertex pairs within r, as in greedy_vertex_walk. Such a walk gives d_F(P, Q) <= r;
# false says nothing. lp and lq are the prefix lengths of the curves.
function greedy_ball_walk(p, q, lp, lq, r::Float64)::Bool
    n, m = length(p), length(q)
    r_sq = r * r
    # The margins cover the rounding of the arc lengths.
    r_in = r * (1.0 - 1.0e-9) - 1.0e-11 * (lp[end] + lq[end])
    i = j = 1
    step = max(n, m)
    while (i < n) || (j < m)
        if step <= 1
            i, j = vertex_step(p, q, i, j, r_sq)
            i == 0 && return false
            step = 2
            continue
        end
        i2, j2 = ball_step_target(n, m, i, j, step)
        if balls_within(piece_ball(p, lp, i, i2), piece_ball(q, lq, j, j2), r_in)
            i, j = i2, j2
            step = cld(3 * step, 2)
        else
            step ÷= 2
        end
    end
    return true
end

# The target of a step of `step` vertices from (i, j): on the curve that is not done,
# or on both, at rates proportional to what remains of them.
function ball_step_target(n::Int, m::Int, i::Int, j::Int, step::Int)
    i == n && return n, min(j + step, m)
    j == m && return min(i + step, n), m
    return min(i + step, n), min(j + max(1, step * (m - j) ÷ (n - i)), m)
end

# The closest of the next vertex pairs from (i, j), if within the radius; (0, 0)
# otherwise.
function vertex_step(p, q, i::Int, j::Int, r_sq::Float64)
    n, m = length(p), length(q)
    d_both = (i < n && j < m) ? fp_dist_sq(p[i+1], q[j+1]) : Inf
    d_p = i < n ? fp_dist_sq(p[i+1], q[j]) : Inf
    d_q = j < m ? fp_dist_sq(p[i], q[j+1]) : Inf
    d = min(d_both, d_p, d_q)
    d <= r_sq || return 0, 0
    d == d_both && return i + 1, j + 1
    return d == d_p ? (i + 1, j) : (i, j + 1)
end

# The squared distance from x to the box of node v of t.
function node_gap_sq(t::AbstractBoxTree, v::Int, x)::Float64
    node_box!(t, v)
    gap = 0.0
    B = t.box
    d = size(B, 1) >> 1
    # In bounds: v is a node of t, and x has d coordinates.
    @inbounds for c in 1:d
        g = max(0.0, B[c, v] - x[c], x[c] - B[d+c, v])
        gap += g * g
    end
    return gap
end

# Whether some point of the curve q (with tree t) is within the radius of x: a descent of
# the tree, nearer child first, that skips the nodes whose box is farther than the radius
# (with a margin), and decides on single segments by their free interval, as the sweep
# does.
function near_curve(x, q, t::PieceTree, r_sq::Float64, v::Int = 1)::Bool
    rng = node_range(t, v)
    isempty(rng) && return false
    node_gap_sq(t, v, x) > r_sq * (1.0 + 1.0e-9) && return false
    is_bottom(t, v) &&
        return any(c -> !free_is_empty(free_interval(x, q[c], q[c+1], r_sq)), rng)
    a, b = 2v, 2v + 1
    node_gap_sq(t, b, x) < node_gap_sq(t, a, x) && ((a, b) = (b, a))
    return near_curve(x, q, t, r_sq, a) || near_curve(x, q, t, r_sq, b)
end

# Whether vertex i of n is within the radius of the vertex of q at the same relative
# position, or of one of its neighbors (a quick sufficient test for near_curve).
function near_guess(x, q, i::Int, n::Int, r_sq::Float64)::Bool
    m = length(q)
    k = 1 + round(Int, (i - 1) * (m - 1) / (n - 1))
    return any(t -> fp_dist_sq(x, q[t]) <= r_sq, max(1, k-1):min(m, k+1))
end

# Whether some vertex of p, among a few spread along it (at offsets growing by a factor
# 1.5 from either end), has no point of q within the radius; then d_F(P, Q) > r.
function far_vertex(p, q, tq::PieceTree, r_sq::Float64)::Bool
    n = length(p)
    step = 1
    while step < n
        for i in (1 + step, n - step)
            (1 < i < n) &&
                !near_guess(p[i], q, i, n, r_sq) &&
                !near_curve(p[i], q, tq, r_sq) &&
                return true
        end
        step = max(step + 1, (3 * step) ÷ 2)
    end
    return false
end
