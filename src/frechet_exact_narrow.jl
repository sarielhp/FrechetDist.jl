# The exact distance from the narrow interval left by a bisection with the fast decision
# procedure: the critical values in it, as list_critical! finds them there (bit for bit),
# computed from the coordinates segment by segment, with the kinetic sweep only on the
# segments where two points can change order; and the search of exact_without_tables,
# which decides with the fast decision procedure instead of the tables.

# The critical values in [lo, hi] of the segments of p against the vertices of q, as
# list_critical! finds them, for a narrow [lo, hi] (that of a bisection): from the
# coordinates, segment by segment (narrow_segment!), and with the kinetic
# sweep only on the segments that are not quiet (narrow_segment!).
function narrow_critical!(out, p, q, lo::Float64, hi::Float64)
    tab = EdgeTable(length(q))
    sw = MonoSweep(tab)
    reach = hi + SWEEP_COINCIDENT * hi
    # Squared radii, and the windows (relative margin 1e-12) outside which a value is
    # certainly not in [lo, hi], and a point certainly not entering or leaving the
    # segment in [lo, reach].
    w = (;
        lo2 = lo * lo,
        r2 = reach * reach,
        win = (lo * lo * (1 - 1e-12), hi * hi * (1 + 1e-12)),
        span = (lo * lo * (1 - 1e-12), reach * reach * (1 + 1e-12)),
        lo,
        hi,
    )
    pos = Float64[]
    for i in 1:(length(p)-1)
        narrow_segment!(out, tab, pos, (p, q, i), w) || mono_values!(out, sw, lo, hi)
    end
    return out
end

# Segment i of p against the vertices of q: fills `tab` with the values edge_table!
# copies from the tables, computed by the same functions (so the same, bit for bit);
# pushes the vertex-vertex and vertex-segment values in [lo, hi]; and returns whether the
# segment is quiet: certainly without a monotonicity value in [lo, reach], so that the
# sweep would record nothing. It is quiet when every point (pt_span) is on the segment
# throughout [lo, reach] or off it throughout (a point near the ends of the window
# counts as an event), and adjacent points on it at lo are more than SWEEP_POS_SLACK +
# 2 delta apart, delta a bound on how far any of them moves from lo to reach. Then they
# are in the same order at reach, still more than SWEEP_POS_SLACK apart; and two points
# cross at most once (their distance along the segment is monotone in the radius), so
# none crosses. `pos` is a buffer.
function narrow_segment!(out, tab::EdgeTable, pos::Vector{Float64}, (p, q, i), w)
    a, b = p[i], p[i+1]
    tab.L2 = exact_dist_sq(a, b)
    tab.L = sqrt(tab.L2)
    # Three passes: the first, alone, is the tightest loop.
    for k in eachindex(q)
        tab.t0[k], tab.h2[k] = line_params(q[k], a, b)
        tab.da2[k], tab.db2[k] = exact_dist_sq(q[k], a), exact_dist_sq(q[k], b)
    end
    last = i == length(p) - 1
    for k in eachindex(q)
        narrow_values!(out, (tab.t0[k], tab.h2[k], tab.da2[k], tab.db2[k]), last, w)
    end
    empty!(pos)
    delta = 0.0
    for k in eachindex(q)
        delta = quiet_points!(
            pos,
            delta,
            (tab.t0[k], tab.h2[k], tab.da2[k], tab.db2[k]),
            tab.L,
            w,
        )
        delta < 0 && return false
    end
    gap = SWEEP_POS_SLACK + 2 * delta
    sort!(pos)
    for j in 2:length(pos)
        pos[j] - pos[j-1] > gap || return false
    end
    return true
end

# The values of edge_point_values! for one vertex, in [lo, hi]; the square root is taken
# only of the squares in the window w.win around [lo^2, hi^2].
@inline function narrow_values!(out, (t0, h2, da2, db2), last::Bool, w)
    near(x2) = w.win[1] <= x2 <= w.win[2]
    vs2 = (0 <= t0 <= 1) ? h2 : min(da2, db2)
    near(vs2) && in_range_push!(out, vertex_segment_value(t0, h2, da2, db2), w.lo, w.hi)
    near(da2) && in_range_push!(out, sqrt(da2), w.lo, w.hi)
    last && near(db2) && in_range_push!(out, sqrt(db2), w.lo, w.hi)
    return nothing
end

# The two points of one vertex for the quiet test of narrow_segment!: pushes the
# positions at lo of those on the segment throughout, and returns the bound delta, or -1
# when one of them enters or leaves the segment in (or near) the window.
@inline function quiet_points!(pos, delta::Float64, (t0, h2, da2, db2), L::Float64, w)
    h2 > w.span[2] && return delta  # its sphere does not reach the segment's line
    s = sqrt(max(w.lo2 - h2, 0.0))
    for sg in (-1.0, 1.0)
        ends = quiet_span(sg, t0, h2, da2, db2)
        ends === nothing && continue
        (ends[1] > w.span[2] || ends[2] < w.span[1]) && continue  # off throughout
        (ends[1] < w.span[1] && ends[2] > w.span[2] && s > 0) || return -1.0
        # sqrt(r2 - h2) - s = (r2 - lo2) / (sqrt(r2 - h2) + s) <= (r2 - lo2) / 2 s.
        delta = max(delta, (w.r2 - w.lo2) / (2 * L * s) * (1 + 1e-3))
        push!(pos, t0 + sg * s / L)
    end
    return delta
end

# The squared radii at which the point behind (sg < 0) or ahead of the closest point of
# a vertex enters and leaves the segment, as pt_span (in squares); nothing if never.
@inline function quiet_span(
    sg::Float64,
    t0::Float64,
    h2::Float64,
    da2::Float64,
    db2::Float64,
)
    if sg < 0
        t0 < 0 && return nothing
        return (t0 > 1 ? db2 : h2, da2)
    end
    t0 > 1 && return nothing
    return (t0 < 0 ? da2 : h2, db2)
end

########################################################################
# The bisection, and the exact distance without tables.

# The curves prepared for the fast decision procedure, without their widths (a pass over
# the curves that a bisection rarely needs).
exact_decide_curves(P, Q) = (
    DecideCurve(P, default_tree(P.pnts); width = false),
    DecideCurve(Q, default_tree(Q.pnts); width = false),
)

# Whether the bisection replaces the sample: long curves in low dimension
# (EXACT_BISECT_FACTOR).
function exact_bisects(p, q)
    n, m = length(p), length(q)
    (n >= 2 && m >= 2) || return false
    return max(n, m)^2 >= EXACT_BISECT_FACTOR * (n + m) * length(p[1])
end

# [lo, hi] narrowed by the fast decision procedure to [a, b], b - a at most
# EXACT_BISECT_RTOL a (or no Float64 between them).
function exact_bisect(A::DecideCurve, B::DecideCurve, lo::Float64, hi::Float64)
    a, b = lo, hi
    while b - a > EXACT_BISECT_RTOL * a
        # Halve the ratio while it is large, then the difference.
        mid = (a > 0.0 && b > 4.0 * a) ? sqrt(a * b) : (a + b) / 2
        (a < mid < b) || break
        frechet_c_decide_exact(A, B, mid) ? (b = mid) : (a = mid)
    end
    return a, b
end

# The critical values in (a, top), sorted and distinct, and top: b plus 16 times the
# slack in ulps (the accepted critical value can lie some ulps above b). The listing
# starts a little below a, as exact_search's does.
function narrow_candidates(P, Q, a::Float64, b::Float64, slack::Int)
    top = nextfloat(b, 16 * slack)
    from = prevfloat(a, 4 * slack)
    out = narrow_critical!(Float64[], P.pnts, Q.pnts, from, top)
    narrow_critical!(out, Q.pnts, P.pnts, from, top)
    return filter!(v -> a < v < top, unique!(sort!(out))), top
end

# The exact distance with the fast decision procedure in place of the tables: O(n + m)
# memory, and decisions that search only the reachable part of the free space. The same
# search as exact_from_bisection: the decider (with the same slack) must reject a and
# accept b; the result is the smallest accepted critical value in (a, top), and nothing
# well below it may be accepted. Returns nothing when a check fails; the caller then
# computes with the tables.
function exact_without_tables(P, Q)
    p, q = P.pnts, Q.pnts
    A, B = exact_decide_curves(P, Q)
    slack = EXACT_SLACK_ULPS
    accepts(c) = frechet_c_decide_exact(A, B, nextfloat(c, slack))
    lo = sqrt(max(exact_dist_sq(p[1], q[1]), exact_dist_sq(p[end], q[end])))
    accepts(lo) && return lo  # a lower bound, and a critical value
    a, b = exact_bisect(A, B, lo, sqrt(box_far_sq(A, B)) * (1.0 + 1.0e-12))
    (accepts(a) || !accepts(b)) && return nothing
    S, _ = narrow_candidates(P, Q, a, b, slack)
    i, j = 0, length(S) + 1
    while j - i > 1
        mid = (i + j) >>> 1
        accepts(S[mid]) ? (j = mid) : (i = mid)
    end
    j > length(S) && return nothing
    x = S[j]
    frechet_c_decide_exact(A, B, prevfloat(x, slack)) && return nothing
    return x
end
