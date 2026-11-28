# The continuous Fréchet distance to a given precision, by bisection with the decision
# procedure (frechet_c_decide_exact), as FrechetLight::calcDistance of Bringmann,
# Künnemann and Nusser does with theirs: most of the decisions are settled by the filters
# and the boundary rules of the free-space search, and only the last ones, close to the
# distance, search the free space at all.

"""
    frechet_c_compute_bisect(P, Q; rtol = 1e-10, atol = 0.0) -> Float64

The continuous Fréchet distance between the curves `P` and `Q`, to within
`max(atol, rtol * d) / 2` (`d` the distance): a bisection with the exact decision
procedure (`frechet_c_decide_exact`) on prepared curves (`DecideCurve`), until the
interval known to hold the distance is at most `max(atol, rtol * lo)` long; the result is
its middle.

The interval starts at the distance of the end points, a lower bound, and at the largest
distance between the bounding boxes of the curves, an upper bound. The first decision is
at the lower bound plus the tolerance: the distance is often that of the end points (the
curves then end at their farthest pair), and then one decision settles it. The curves
are prepared without their widths (`DecideCurve(...; width = false)`): computing them
costs a pass over the curves that the decisions of a bisection rarely need.

Each decision costs O(n + m) for its filters and, when they do not settle it, a search
of the free space by boxes that is usually much faster than the O(n m) of the cell by
cell sweep; about log2((hi - lo) / tolerance) decisions are made. With `rtol` down to
about 1e-14 the result is the distance up to the rounding of the decision procedure;
`frechet_c_compute_exact` gives the critical value itself.

`atol = 1e-10, rtol = 0` is the stopping rule of FrechetLight::calcDistance
(the C++ code of Bringmann et al.).
"""
function frechet_c_compute_bisect(
    P::GPntSeq{N,T,PT},
    Q::GPntSeq{N,T,PT};
    rtol::Real = 1e-10,
    atol::Real = 0.0,
)::Float64 where {N,T,PT<:AbstractVector{T}}
    who = "frechet_c_compute_bisect"
    check_curves(P, Q, who; min_vertices = 1)
    rt = check_tolerance(rtol, who)
    at = check_tolerance(atol, who)
    (rt > 0.0 || at > 0.0) ||
        throw(ArgumentError("$who: one of rtol and atol must be positive"))
    p, q = P.pnts, Q.pnts
    # A curve of one vertex: the farthest vertex of the other curve.
    length(p) == 1 && return sqrt(maximum(y -> fp_dist_sq(p[1], y), q))
    length(q) == 1 && return sqrt(maximum(x -> fp_dist_sq(x, q[1]), p))
    A = DecideCurve(P, default_tree(p); width = false)
    B = DecideCurve(Q, default_tree(q); width = false)
    lo = sqrt(max(fp_dist_sq(p[1], q[1]), fp_dist_sq(p[end], q[end])))
    hi = sqrt(box_far_sq(A, B)) * (1.0 + 1.0e-12)
    tol = max(at, rt * lo)
    hi - lo <= tol && return (lo + hi) / 2
    probe = lo + tol
    frechet_c_decide_exact(A, B, probe) && return (lo + probe) / 2
    lo = probe
    while hi - lo > max(at, rt * lo)
        # Halve the ratio while it is large, then the difference.
        mid = (lo > 0.0 && hi > 4.0 * lo) ? sqrt(lo * hi) : (lo + hi) / 2
        (lo < mid < hi) || break  # no Float64 between them
        frechet_c_decide_exact(A, B, mid) ? (hi = mid) : (lo = mid)
    end
    return (lo + hi) / 2
end
