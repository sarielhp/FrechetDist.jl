# The pieces of the free-space search (frechet_decide_blocks.jl) as sausages instead of
# bounding boxes: a piece of the curve, from vertex a to vertex b, lies in the Minkowski
# sum of the segment ab and the ball of radius rho around the origin, rho the largest
# distance from a vertex of the piece to ab (the sum is convex, so it contains the
# segments between the vertices too). Unlike a bounding box, a sausage does not grow
# with the dimension, and its tests need only distances between points:
#
# - the distance from a point to a segment, from the squared distances from the point to
#   its ends (the inner product by polarization);
# - two pieces are within r of each other (a free box) when the farthest end points of
#   their segments are within r - rho1 - rho2 (the distance between points of two
#   segments is convex, so it is largest at end points).
#
# A node also stores phi >= rho, an upper bound on the Fréchet distance between its
# piece and its segment, for greedy_sausage_walk (corner to corner, where rho says
# nothing: a piece going back and forth over its segment has rho = 0).
#
# The tree has the nodes of BoxTree (same numbering and ranges); a node stores rho and
# phi (its segment joins the end vertices of its piece), computed the first time the
# node is read: by a scan of its piece at a bottom node, from its children above
# (fill_radii!; reading the root computes the whole tree, in O(n d) time). The radii are
# computed from coordinate differences and enlarged by radius_margin, which covers their
# rounding; the distances computed by polarization (node_gap_sq) are corrected by
# SAUSAGE_MARGIN times the largest squared distance involved; and the tests have the
# relative margin 1e-9 of the box tests: they decide a box only when the cell by cell
# computation would agree. Reading a node writes into the tree, so a tree must not be
# used by several threads at once.
struct SausageTree{V} <: PieceTree
    p::V
    # rad[1, v] = rho and rad[2, v] = phi of node v: NaN until computed, -1 for the nodes
    # past the end of the curve.
    rad::Matrix{Float64}
    N::Int
    ncell::Int
    lens::Vector{Float64}  # empty until tree_lens is first called
end

function SausageTree(p)
    ncell = length(p) - 1
    N = nextpow(2, max(cld(ncell, DECIDE_BLOCK), 1))
    return SausageTree(p, fill(NaN, 2, 2N), N, ncell, Float64[])
end

tree_kind(::SausageTree) = SausageTree

"Relative margin, to the largest squared distance involved, of the sausage distances."
const SAUSAGE_MARGIN = 1.0e-11
tree_lens(t::SausageTree) = lazy_lens!(t.lens, t.p)

@inline function node_box!(t::SausageTree, v::Int)
    # In bounds: v is a node of t, so v <= 2N.
    isnan(@inbounds t.rad[1, v]) && fill_radii!(t, v)
    return nothing
end

# The squared distance from a point to a segment, from the squared distances A and B
# from the point to its ends and the squared length L of the segment.
function seg_dist_sq(A::Float64, B::Float64, L::Float64)::Float64
    L > 0.0 || return A
    s = clamp((A + L - B) / (2.0 * L), 0.0, 1.0)
    return max(0.0, A - s * (A + L - B) + s * s * L)
end


# The segment ab of a piece, for the distances from the vertices of the piece to it.
struct SegFrame{P}
    a::P
    b::P
    L::Float64  # |ab|^2
    lab::Float64  # |ab|
    d::Int
end

function SegFrame(p, i::Int, k::Int)
    a, b = p[i], p[k]
    L = fp_dist_sq(a, b)
    return SegFrame(a, b, L, sqrt(L), box_dim(p))
end

# The parameter of the projection of x on the line ab, clamped to [0, 1], and |xa|^2.
@inline function proj_param(x, f::SegFrame)
    a, b = f.a, f.b
    t = s = 0.0
    # In bounds: x, a and b have d coordinates.
    @inbounds for e in 1:f.d
        xa = x[e] - a[e]
        t += xa * (b[e] - a[e])
        s += xa * xa
    end
    return (f.L > 0.0 ? clamp(t / f.L, 0.0, 1.0) : 0.0), s
end

# The distance from x to the point a + tau (b - a), from the coordinate differences
# ((x - a) - tau (b - a), exact near 0, unlike polarization).
@inline function resid(x, f::SegFrame, tau::Float64)::Float64
    a, b = f.a, f.b
    s = 0.0
    # In bounds: x, a and b have d coordinates.
    @inbounds for e in 1:f.d
        y = (x[e] - a[e]) - tau * (b[e] - a[e])
        s += y * y
    end
    return sqrt(s)
end

# A bound on the rounding error of resid, for |xa| and |ab| at most len: the coordinates
# of the difference vector are within about 3u (|xa| + |ab|) <= 6u len of the exact
# ones in norm, and the sum of d squares and the root add (d + 3) u times the result, at
# most 2 len (u = eps/2). So the error is below (d + 9) u len; twice that is the margin.
radius_margin(d::Int, len::Float64) = (d + 9) * eps(Float64) * len

# The state of a scan of the vertices of a piece against its segment ab:
# (tau, phi, rho, low). Each vertex x is matched to the point of ab at the largest
# parameter tau of the projections so far (a monotone matching, along which the distance
# is convex on each edge, so largest at a vertex): phi, the largest distance in this
# matching, bounds the Fréchet distance of the piece to ab from above. rho is the
# largest distance from a vertex to ab (from above); low, the same from below.
const SCAN_START = (0.0, 0.0, 0.0, 0.0)

@inline function scan_point(x, f::SegFrame, st::NTuple{4,Float64})
    tau, phi, rho, low = st
    tc, xa_sq = proj_param(x, f)
    margin = radius_margin(f.d, max(sqrt(xa_sq), f.lab))
    near = resid(x, f, tc)
    tc > tau && (tau = tc)
    far = tau == tc ? near : resid(x, f, tau)
    return (tau, max(phi, far + margin), max(rho, near + margin), max(low, near - margin))
end

# The scan of the vertices p[c] for c in cs.
function scan_points(p, cs::UnitRange{Int}, f::SegFrame, st::NTuple{4,Float64})
    for c in cs
        st = scan_point(p[c], f, st)
    end
    return st
end

# rho and phi of node v. At a bottom node, by a scan of its piece. Above, from the
# children, split at the vertex m: rho = max(rho_left, rho_right) + d(m, ab) (the
# sausage of the left child, around am, is within d(m, ab) + rho_left of ab, as the
# distance to ab is convex; the same on the right), and the same for phi (the matchings
# of the children concatenate: the piece is within max(phi_left, phi_right) of the
# polygon amb in Fréchet distance, and d_F(amb, ab) = d(m, ab): match m to its nearest
# point of ab). So a node above the bottom costs O(d), not a scan of its piece: the
# sums are loose by a constant factor on rough curves, but in tests the search was
# faster with them than with exact radii by rescans (in the plane, by up to 1.3 times).
@noinline function fill_radii!(t::SausageTree, v::Int)
    rng = node_range(t, v)
    if isempty(rng)
        t.rad[1, v] = t.rad[2, v] = -1.0
        return nothing
    end
    i, k = first(rng), last(rng) + 1
    f = SegFrame(t.p, i, k)
    if is_bottom(t, v)
        _, phi, rho, _ = scan_points(t.p, (i+1):(k-1), f, SCAN_START)
        t.rad[1, v], t.rad[2, v] = rho, phi
        return nothing
    end
    node_box!(t, 2v)
    node_box!(t, 2v + 1)
    if t.rad[1, 2v+1] < 0.0  # the right child is past the end: v has the left's piece
        t.rad[1, v], t.rad[2, v] = t.rad[1, 2v], t.rad[2, 2v]
        return nothing
    end
    _, _, dm, _ = scan_point(t.p[first(node_range(t, 2v + 1))], f, SCAN_START)
    t.rad[1, v] = max(t.rad[1, 2v], t.rad[1, 2v+1]) + dm
    t.rad[2, v] = max(t.rad[2, 2v], t.rad[2, 2v+1]) + dm
    return nothing
end

# The end vertices of the piece of node v (whose range is not empty).
node_ends(t::SausageTree, v::Int) =
    (rng = node_range(t, v); (t.p[first(rng)], t.p[last(rng)+1]))

# Whether every point of the piece of node u is within the radius of every point of the
# piece of node v: the distance between points of two segments is convex, so it is
# largest at end points, and the sausages add their radii.
function pieces_within(
    tp::SausageTree,
    tq::SausageTree,
    u::Int,
    v::Int,
    r_sq::Float64,
)::Bool
    node_box!(tp, u)
    node_box!(tq, v)
    a1, b1 = node_ends(tp, u)
    a2, b2 = node_ends(tq, v)
    ends_sq =
        max(fp_dist_sq(a1, a2), fp_dist_sq(a1, b2), fp_dist_sq(b1, a2), fp_dist_sq(b1, b2))
    far = sqrt(ends_sq) + tp.rad[1, u] + tq.rad[1, v]
    return far * far * (1.0 + 1.0e-9) <= r_sq
end

# A lower bound on the squared distance from x to the piece of node v.
function node_gap_sq(t::SausageTree, v::Int, x)::Float64
    node_box!(t, v)
    t.rad[1, v] < 0.0 && return Inf
    a, b = node_ends(t, v)
    A, B, L = fp_dist_sq(x, a), fp_dist_sq(x, b), fp_dist_sq(a, b)
    g = sqrt(max(0.0, seg_dist_sq(A, B, L) - SAUSAGE_MARGIN * max(A, B, L))) - t.rad[1, v]
    return g > 0.0 ? g * g : 0.0
end

# An upper bound on the Fréchet distance between the piece p[i:k] and the segment from
# p[i] to p[k], by the matching of scan_point. By one scan.
segment_frechet(p, i::Int, k::Int)::Float64 =
    scan_points(p, (i+1):(k-1), SegFrame(p, i, k), SCAN_START)[2]

# Bounds on segment_frechet(p, i, k) from the tree, in O((DECIDE_BLOCK + log n) d) time,
# as (upper, lower). p[i:k] is the vertices up to the first bottom node boundary, a run
# of whole nodes (the largest aligned ones), and the vertices after the last boundary.
# Upper: the Fréchet distance of the piece to the polygon through those vertices and the
# ends of the nodes is at most the largest phi of the nodes (the matchings concatenate),
# and that of the polygon to the segment is bounded by the matching of scan_point; the
# two add up. Lower: every vertex of the polygon is a vertex of the piece, and the scan
# bounds the Fréchet distance from above, so from below by the distance of a vertex to
# the segment. (Short pieces are scanned: both bounds are the scan.)
function frechet_bounds(t::SausageTree, i::Int, k::Int)
    B, p = DECIDE_BLOCK, t.p
    s0, s1 = cld(i - 1, B), (k - 1) ÷ B  # the bottom nodes s0 .. s1-1 lie in p[i:k]
    s1 - s0 < 2 && (u = segment_frechet(p, i, k); return (u, u))
    f = SegFrame(p, i, k)
    st = scan_points(p, (i+1):(s0*B+1), f, SCAN_START)
    phi = 0.0
    s = s0
    while s < s1
        size = s == 0 ? t.N : 1 << trailing_zeros(s)  # the largest aligned node at s
        while size > s1 - s
            size >>= 1
        end
        v = (t.N + s) >> trailing_zeros(size)
        node_box!(t, v)
        phi = max(phi, t.rad[2, v])
        s += size
        st = scan_point(p[s*B+1], f, st)
    end
    st = scan_points(p, (s1*B+2):(k-1), f, st)
    return (phi + st[2], st[4])
end

# Where the walk gets the Fréchet bounds of pieces: a SausageTree, or a scan of the
# vertices (ScanPieces) for the other trees.
struct ScanPieces{V}
    p::V
end
walk_pieces(p, ::Nothing) = ScanPieces(p)
walk_pieces(p, ::PieceTree) = ScanPieces(p)
walk_pieces(_, t::SausageTree) = t
frechet_bounds(s::ScanPieces, i::Int, k::Int) = (u = segment_frechet(s.p, i, k); (u, u))
walk_points(s::ScanPieces) = s.p
walk_points(t::SausageTree) = t.p

# Whether d_F(P[is], Q[js]) <= r follows from the sausages of the two pieces: the
# Fréchet distance between two segments is the larger distance between their ends, and
# each piece is within segment_frechet of its segment. By the bounds of frechet_bounds
# first; only when they do not decide, by the scans.
function sausages_within(P, is::UnitRange{Int}, Q, js::UnitRange{Int}, r::Float64)::Bool
    p, q = walk_points(P), walk_points(Q)
    ends = max(fp_dist(p[first(is)], q[first(js)]), fp_dist(p[last(is)], q[last(js)]))
    ends <= r || return false
    up_p, low_p = frechet_bounds(P, first(is), last(is))
    ends + low_p <= r || return false
    up_q, low_q = frechet_bounds(Q, first(js), last(js))
    ends + up_p + up_q <= r && return true
    ends + low_p + low_q <= r || return false
    thick = min(up_p, segment_frechet(p, first(is), last(is)))
    ends + thick + low_q <= r || return false
    return ends + thick + min(up_q, segment_frechet(q, first(js), last(js))) <= r
end


"""Smallest dimension in which the decision tries greedy_sausage_walk: in lower dimensions its
scans cost more than they settle (the boxes and balls describe the pieces well)."""
const SAUSAGE_WALK_DIM = 6

# The walk of greedy_ball_walk, by sausages instead of balls: from (i, j) to (i', j')
# when sausages_within shows d_F(p[i:i'], q[j:j']) <= r. With SausageTrees tp and tq a
# step costs O((DECIDE_BLOCK + log n) d) (their phi); otherwise a scan of the two pieces,
# so a walk costs O((n + m) d) up to the failed attempts. It is tight on pieces that are
# nearly straight, in any dimension (two close segments, refined into many vertices, are
# settled by the first step). True gives d_F(P, Q) <= r; false says nothing.
function greedy_sausage_walk(p, q, r::Float64, tp = nothing, tq = nothing)::Bool
    P, Q = walk_pieces(p, tp), walk_pieces(q, tq)
    n, m = length(p), length(q)
    r_sq = r * r
    r_in = r * (1.0 - 1.0e-9)  # covers the rounding of the sums of distances
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
        if sausages_within(P, i:i2, Q, j:j2, r_in)
            i, j = i2, j2
            step = cld(3 * step, 2)
        else
            step ÷= 2
        end
    end
    return true
end

"""Dimension from which the decision describes the pieces of the curves by sausages
(SausageTree) instead of bounding boxes (LazyBoxTree), by default: below it boxes are
tight and cheap; above it they grow with the dimension, and sausages do not."""
const SAUSAGE_TREE_DIM = 64

"The kind of tree the decision uses for the vertices p, by default (see SAUSAGE_TREE_DIM)."
default_tree(p) =
    (!isempty(p) && box_dim(p) >= SAUSAGE_TREE_DIM) ? SausageTree : LazyBoxTree
