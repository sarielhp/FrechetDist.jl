# A BoxTree (frechet_decide_blocks.jl) built on demand: the same implicit tree (bottom
# nodes of DECIDE_BLOCK segments, children 2v and 2v + 1, parent v >> 1) and the same
# layout (column v of box the bounding box of node v), but a box is computed the first
# time it is read, and the prefix lengths the first time they are needed:
#
# - the first read builds the top of the tree, the nodes of at least LAZY_PIECE
#   segments, by one scan of the curve (the boxes of its pieces of LAZY_PIECE segments,
#   then their unions);
# - the first read of a smaller node builds its whole subtree, by one scan of its
#   bottom nodes and their unions.
#
# So each vertex is read at most about twice, and below the top only the pieces the
# search reaches are built. A decision settled by the filters reads nothing; one settled
# near the root reads only the top. A box not built yet has a NaN first coordinate, so
# the check reads the cache line the box is read from anyway.
#
# The boxes are minima and maxima of coordinates, so they are exactly those of BoxTree,
# and so are the decisions. Reading a node writes into the tree, so a tree (and a
# DecideCurve holding it) must not be used by several threads at once.
struct LazyBoxTree{V} <: AbstractBoxTree
    p::V
    box::Matrix{Float64}  # 2d x 2N, as in BoxTree; NaN in row 1 until built
    N::Int
    ncell::Int
    lens::Vector{Float64}  # empty until tree_lens is first called
end

function LazyBoxTree(p)
    ncell = length(p) - 1
    N = nextpow(2, max(cld(ncell, DECIDE_BLOCK), 1))
    box = fill(NaN, 2 * box_dim(p), 2N)
    return LazyBoxTree(p, box, N, ncell, Float64[])
end

tree_kind(::LazyBoxTree) = LazyBoxTree

@inline function node_box!(t::LazyBoxTree, v::Int)
    # In bounds: v is a node of t, so v < 2N.
    isnan(@inbounds t.box[1, v]) && fill_node_box!(t, v)
    return nothing
end

"Size (in segments) of the pieces at the bottom of the top part of a LazyBoxTree."
const LAZY_PIECE = 128

@noinline function fill_node_box!(t::LazyBoxTree, v::Int)
    if isnan(t.box[1, 1])
        piece = min(t.N, prevpow(2, max(1, LAZY_PIECE ÷ DECIDE_BLOCK)))  # bottom nodes per piece
        fill_levels!(t, 1, trailing_zeros(t.N) - trailing_zeros(piece))
        isnan(t.box[1, v]) || return nothing
    end
    fill_levels!(t, v, trailing_zeros(node_size(t, v) ÷ DECIDE_BLOCK))
    return nothing
end

# The boxes of node v and of its descendants down to k levels below it: those k levels
# below from their vertices, the others from their children.
function fill_levels!(t::LazyBoxTree, v::Int, k::Int)
    for w in (v<<k):((v<<k)+(1<<k)-1)
        scan_box!(t.box, t.p, node_range(t, w), w)
    end
    d = box_dim(t.p)
    for h in (k-1):-1:0, w in (v<<h):((v<<h)+(1<<h)-1)
        union_box!(t.box, w, d)
    end
    return nothing
end

function tree_lens(t::LazyBoxTree)
    lens, p = t.lens, t.p
    isempty(lens) || return lens
    resize!(lens, t.ncell + 1)
    lens[1] = 0.0
    for c in 1:(t.ncell)
        lens[c+1] = lens[c] + fp_dist(p[c], p[c+1])
    end
    return lens
end
