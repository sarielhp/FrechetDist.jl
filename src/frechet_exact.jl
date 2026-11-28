# The exact continuous Fréchet distance, by the randomized search over critical values
# of Har-Peled and Raichel, "The Fréchet distance revisited and extended" (Section 6,
# "Computing optimal Fréchet paths for DAG complexes", for the case of two curves).
#
# The distance is one of the critical values of Alt and Godau (1995), where the answer
# of the decision procedure can change:
#   - vertex-vertex: |p_i - q_j|;
#   - vertex-segment: the distance of a vertex of one curve to a segment of the other;
#   - monotonicity: for a segment e of one curve and vertices q_k, q_l (k < l) of the
#     other, the radius at which their spheres meet on e (the point of e equidistant
#     from both, if it lies on e).
# There are O(n^3) of them. The algorithm binary searches a random sample of 4 N^2 of
# them with the decider, which leaves an interval [a, b] holding the distance and, with
# high probability, only O(N log N) critical values. Those are listed explicitly (the
# monotonicity values by a kinetic sweep along each segment), sorted, and searched
# again.
#
# Everything taken from the coordinates is computed once, into tables (`ExactTables`,
# O(n m d) time): for every vertex and segment of the other curve, the closest point on
# the line and the squared distance to it (`line_params`), and the squared vertex-vertex
# distances (`exact_dist_sq`). Sampling, listing and every decision of the search then
# cost O(1) per pair of vertex and segment, instead of O(d). The computations are plain
# loops (no `@simd` or `@fastmath`, which let the compiler reorder the sums differently
# where they are inlined), and the critical values come from the tables only, so the
# same value sampled and listed is the same Float64.

# A candidate radius c counts as accepted when the decider accepts c plus this many
# ulps. The decider and the critical values read the same tables, but combine them
# differently, so the radius where the decider flips lies some ulps on either side of
# the computed critical value. That error does not grow with the dimension: no sum
# over coordinates is involved after the tables (measured: up to 16 ulps on real
# curves in R^963, 1 in R^19200; deciding from the coordinates instead, the offsets
# grew to 162 ulps in R^19200).
const EXACT_SLACK_ULPS = 64

# How often the result failed the minimality check and the critical values were listed
# by brute force instead (a missed value in the sweep). Read by the tests.
const EXACT_FALLBACKS = Threads.Atomic{Int}(0)

"Squared distance of two points, summed in order (see the note at the top)."
function exact_dist_sq(x::AbstractVector{T}, y::AbstractVector{T}) where {T}
    n = vec_length(x, y)
    s = zero(T)
    # In bounds: x and y have n coordinates, indexed from 1 (vec_length).
    @inbounds for k in 1:n
        s += (x[k] - y[k])^2
    end
    return s
end

"""
    line_params(q, a, b) -> (t0, h2)

For the vertex `q` and the segment `a b`: the parameter `t0` of the closest point on the
line through `a` and `b` (`a + t0 (b - a)`), and the squared distance `h2` to it (from a
second pass, not by cancellation). Assumes `a != b`.
"""
@inline function line_params(
    q::AbstractVector{T},
    a::AbstractVector{T},
    b::AbstractVector{T},
) where {T}
    n = vec_length(q, a, b)
    dd = de = zero(T)
    # In bounds: q, a and b have n coordinates, indexed from 1 (vec_length).
    @inbounds for k in 1:n
        d = b[k] - a[k]
        dd += d * d
        de += d * (q[k] - a[k])
    end
    t0 = de / dd
    h2 = zero(T)
    @inbounds for k in 1:n
        g = (a[k] - q[k]) + t0 * (b[k] - a[k])
        h2 += g * g
    end
    return (t0, h2)
end

"The distance of a vertex to a segment, from its `line_params`."
@inline vertex_segment_value(t0, h2, da2, db2) =
    (0 <= t0 <= 1) ? sqrt(h2) : sqrt(min(da2, db2))

"""
    crossing(L2, t0k, h2k, t0l, h2l) -> (ts, delta)

The monotonicity critical value of the vertices k < l (in this order) against a segment
of squared length `L2`: the parameter `ts` of the point of the segment equidistant from
both, and that distance. `(NaN, NaN)` if there is no such point on the segment (or the
two vertices are equidistant from every point of it).
"""
@inline function crossing(
    L2::Float64,
    t0k::Float64,
    h2k::Float64,
    t0l::Float64,
    h2l::Float64,
)
    ts, d = crossing_line(L2, t0k, h2k, t0l, h2l)
    return (0 <= ts <= 1) ? (ts, d) : (NaN, NaN)
end

# `crossing` on the whole line: the same values, but no test that the point lies on the
# segment (NaN only when the two vertices are equidistant from every point of it).
@inline function crossing_line(
    L2::Float64,
    t0k::Float64,
    h2k::Float64,
    t0l::Float64,
    h2l::Float64,
)
    dt = t0l - t0k
    dt == 0 && return (NaN, NaN)
    ts = (t0k + t0l) / 2 + (h2l - h2k) / (2 * L2 * dt)
    return (ts, sqrt(L2 * (ts - t0k)^2 + h2k))
end

# The segments of one curve against the vertices of the other. The matrices are indexed
# [vertex of the other curve, segment (or vertex) of this one], or, when F (flipped),
# the other way around: ExactTables flips the table of the segments of q, so that both
# are indexed by q first and the decider reads them in memory order.
struct SegTable{F}
    L2::Vector{Float64}   # L2[i]: squared length of segment i
    t0::Matrix{Float64}   # t0, h2 at (k, i): line_params of vertex k of the other curve
    h2::Matrix{Float64}   #   against segment i
    d2::Matrix{Float64}   # d2 at (k, i): squared distance of vertex k of the other curve
end                       #   to vertex i

# Entry (k, i) of the matrix A of s: vertex k of the other curve, segment or vertex i.
@inline at(::SegTable{false}, A::Matrix{Float64}, k::Int, i::Int) = A[k, i]
@inline at(::SegTable{true}, A::Matrix{Float64}, k::Int, i::Int) = A[i, k]
# The entries of A for segment or vertex i, over all vertices of the other curve.
at_all(::SegTable{false}, A::Matrix{Float64}, i::Int) = view(A, :, i)
at_all(::SegTable{true}, A::Matrix{Float64}, i::Int) = view(A, i, :)
other_count(s::SegTable) = length(at_all(s, s.d2, 1))

# The segments of p against the vertices of q; d2 holds the squared distances as
# SegTable{F} indexes them.
function seg_table(p, q, d2::Matrix{Float64}, ::Val{F}) where {F}
    n, m = length(p), length(q)
    L2 = [exact_dist_sq(p[i], p[i+1]) for i in 1:(n-1)]
    t0 = F ? Matrix{Float64}(undef, n - 1, m) : Matrix{Float64}(undef, m, n - 1)
    h2 = similar(t0)
    # The O(n m d) part of the search. Entries are independent, so the threads change
    # nothing in the result.
    Threads.@threads for i in 1:(n-1)
        for k in 1:m
            x, y = line_params(q[k], p[i], p[i+1])
            F ? (t0[i, k] = x; h2[i, k] = y) : (t0[k, i] = x; h2[k, i] = y)
        end
    end
    return SegTable{F}(L2, t0, h2, d2)
end

"""
    exact_tables(p, q) -> ExactTables

Everything `frechet_c_compute_exact` needs from the coordinates of the vertex vectors
`p` and `q`: the segments of each curve against the vertices of the other
(`pq`, `qp`), in O(n m d) time and 40 n m bytes. It is also a free-space source for
`decide_free_space`, at O(1) per free interval.
"""
struct ExactTables
    pq::SegTable{false}  # segments of p, vertices of q
    qp::SegTable{true}   # segments of q, vertices of p
end

function exact_tables(p, q)
    d2 = Matrix{Float64}(undef, length(q), length(p))  # d2[k, i] = |q_k - p_i|^2
    Threads.@threads for i in eachindex(p)
        for k in eachindex(q)
            d2[k, i] = exact_dist_sq(q[k], p[i])
        end
    end
    return ExactTables(seg_table(p, q, d2, Val(false)), seg_table(q, p, d2, Val(true)))
end

# The free interval of vertex k of the other curve on segment i, as free_interval
# computes it from the coordinates, but from the table.
@inline function table_free_interval(s::SegTable, k::Int, i::Int, r_sq::Float64)
    ends = (at(s, s.d2, k, i) <= r_sq, at(s, s.d2, k, i+1) <= r_sq)
    trivial = free_trivial(s.L2[i], ends)
    trivial === nothing || return trivial
    return free_from_line(at(s, s.t0, k, i), at(s, s.h2, k, i), s.L2[i], ends, r_sq)
end

space_size(t::ExactTables) = (size(t.pq.d2, 2), size(t.pq.d2, 1))
free_p_on_q(t::ExactTables, i, j, r_sq) = table_free_interval(t.qp, i, j, r_sq)
free_q_on_p(t::ExactTables, j, i, r_sq) = table_free_interval(t.pq, j, i, r_sq)
within(t::ExactTables, i, j, r_sq) = t.pq.d2[j, i] <= r_sq

# One segment of a SegTable, for the sweep: the values of every vertex of the other
# curve against it.
mutable struct EdgeTable
    L2::Float64
    L::Float64
    t0::Vector{Float64}
    h2::Vector{Float64}
    da2::Vector{Float64}
    db2::Vector{Float64}
end

EdgeTable(m::Int) = EdgeTable(0.0, 0.0, zeros(m), zeros(m), zeros(m), zeros(m))

# Segment i of `s` (O(m): copies).
function edge_table!(tab::EdgeTable, s::SegTable, i::Int)
    tab.L2 = s.L2[i]
    tab.L = sqrt(tab.L2)
    copyto!(tab.t0, at_all(s, s.t0, i))
    copyto!(tab.h2, at_all(s, s.h2, i))
    copyto!(tab.da2, at_all(s, s.d2, i))
    copyto!(tab.db2, at_all(s, s.d2, i + 1))
    return tab
end

@inline in_range_push!(out, x, lo, hi) = (lo <= x <= hi) && push!(out, x)

########################################################################
# Kinetic sweep along one segment (Section 6.2.1): the spheres of radius delta around
# the vertices of the other curve cut the segment in up to two points each; as delta
# grows from lo to hi, two points of different spheres swap their order exactly at a
# monotonicity value of the two vertices. Point 2k-1 is the one of vertex k behind its
# closest point (t <= t0), point 2k the one ahead.

# How far outside [0, 1] (in segment parameter) a crossing may compute and still be
# taken for one at an end of the segment.
const SWEEP_TS_SLACK = 1e-12

# Relative distance within which a crossing counts as falling at an insertion or removal,
# and how close (in segment parameter) the points checked for one must be.
const SWEEP_COINCIDENT = 1e-9
const SWEEP_POS_SLACK = 1e-7

@inline pt_vertex(id::Int) = (id + 1) >> 1
@inline pt_sign(id::Int) = isodd(id) ? -1.0 : 1.0

# Position of point `id` along the segment at radius d, and its speed (for ties).
@inline function pt_key(tab::EdgeTable, id::Int, d::Float64)
    k, sg = pt_vertex(id), pt_sign(id)
    s = sqrt(max(d * d - tab.h2[k], 0.0))
    speed = s == 0 ? sg * Inf : sg * d / (tab.L * s)
    return (tab.t0[k] + sg * s / tab.L, speed)
end

# Where point `id` enters the segment: at the closest point of its sphere, or at an end.
@inline function pt_entry(tab::EdgeTable, id::Int)
    t0 = tab.t0[pt_vertex(id)]
    return t0 < 0 ? 0.0 : (t0 > 1 ? 1.0 : t0)
end

# The radii during which point `id` lies on the segment; empty (Inf, -Inf) if never.
function pt_span(tab::EdgeTable, id::Int)
    k = pt_vertex(id)
    t0 = tab.t0[k]
    h = sqrt(tab.h2[k])
    if pt_sign(id) < 0
        t0 < 0 && return (Inf, -Inf)
        return (t0 > 1 ? sqrt(tab.db2[k]) : h, sqrt(tab.da2[k]))
    end
    t0 > 1 && return (Inf, -Inf)
    return (t0 < 0 ? sqrt(tab.da2[k]) : h, sqrt(tab.db2[k]))
end

# The points on the segment are kept in a balanced search tree (DataStructures'
# SortedSet, a 2-3 tree): O(log m) to insert, delete, or step to a neighbor. Their order
# changes with the radius, so the set holds slots, not points: slots compare by the
# positions of the points in them at the current radius, and two adjacent points swap
# by exchanging slots, which leaves the tree as it is.

# What the order of the slots depends on: the current radius and the point in each slot.
mutable struct SweepClock
    tab::EdgeTable
    item::Vector{Int}                   # item[s]: the point in slot s
    delta::Float64                      # the current radius
    inserting::Int                      # the slot being inserted (0: none) ...
    insert_key::Tuple{Float64,Float64}  # ... compares by this key
end

struct SweepOrdering <: Base.Order.Ordering
    clock::SweepClock
end

@inline function slot_key(c::SweepClock, s::Int)
    s == c.inserting && return c.insert_key
    return pt_key(c.tab, c.item[s], c.delta)
end

# Slots compare by position, then speed; distinct slots never compare equal (the set
# would take them for one element).
function Base.Order.lt(o::SweepOrdering, a::Int, b::Int)
    ka, kb = slot_key(o.clock, a), slot_key(o.clock, b)
    return ka < kb || (ka == kb && a < b)
end

mutable struct MonoSweep
    clock::SweepClock
    set::SortedSet{Int,SweepOrdering}
    slot::Vector{Int}            # slot[id]: the slot of point id, 0 if absent
    token::Vector{SetSemiToken}  # token[s]: where slot s is in the set
    free::Vector{Int}            # unused slots
    heap::BinaryMinHeap{Tuple{Float64,Int,Int}}  # future swaps (radius, left, right)
end

function MonoSweep(tab::EdgeTable)
    cap = 2 * length(tab.t0)
    clock = SweepClock(tab, zeros(Int, cap), 0.0, 0, (0.0, 0.0))
    return MonoSweep(
        clock,
        SortedSet{Int,SweepOrdering}(SweepOrdering(clock)),
        zeros(Int, cap),
        fill(SetSemiToken(0), cap),
        Int[],
        BinaryMinHeap{Tuple{Float64,Int,Int}}(),
    )
end

# The point after (step = advance) or before (step = regress) point id, or 0.
function sweep_neighbor(sw::MonoSweep, id::Int, step::F) where {F}
    st = step((sw.set, sw.token[sw.slot[id]]))
    (st == pastendsemitoken(sw.set) || st == beforestartsemitoken(sw.set)) && return 0
    return sw.clock.item[deref((sw.set, st))]
end

# Puts point id into a free slot of the set, placed by `key` (O(log m)).
function sweep_add!(sw::MonoSweep, id::Int, key::Tuple{Float64,Float64})
    s = pop!(sw.free)
    sw.clock.item[s] = id
    sw.slot[id] = s
    sw.clock.inserting, sw.clock.insert_key = s, key
    sw.token[s] = push_return_semitoken!(sw.set, s)[2]
    sw.clock.inserting = 0
    return nothing
end

# Removes point id, through its token: no comparisons (O(log m)).
function sweep_drop!(sw::MonoSweep, id::Int)
    s = sw.slot[id]
    delete!((sw.set, sw.token[s]))
    sw.slot[id] = 0
    push!(sw.free, s)
    return nothing
end

# Schedules the swap of the adjacent points u (left) and v, if u is to pass v.
function sweep_schedule!(sw::MonoSweep, u::Int, v::Int)
    ku, kv = pt_vertex(u), pt_vertex(v)
    ku == kv && return nothing  # the two points of one sphere never meet again
    tab = sw.clock.tab
    k, l = minmax(ku, kv)
    # Points on the segment can cross only on it, but a crossing at an end (two points
    # entering there together) computes up to rounding outside: SWEEP_TS_SLACK. Only
    # crossings `crossing` puts on the segment are recorded (sweep_swap!).
    ts, d = crossing_line(tab.L2, tab.t0[k], tab.h2[k], tab.t0[l], tab.h2[l])
    (!isnan(d) && -SWEEP_TS_SLACK <= ts <= 1 + SWEEP_TS_SLACK) || return nothing
    # Each point reaches only its own side of its sphere's closest point.
    (
        pt_sign(u) * (ts - tab.t0[ku]) >= -SWEEP_TS_SLACK &&
        pt_sign(v) * (ts - tab.t0[kv]) >= -SWEEP_TS_SLACK
    ) || return nothing
    # The two meet once; the faster one at that radius is behind before and ahead after.
    # Deciding the direction this way, not by comparing d with the current radius,
    # keeps near-simultaneous crossings (radii an ulp apart) in order: a swap due an ulp
    # in the past is still made, next.
    pt_key(tab, u, d)[2] > pt_key(tab, v, d)[2] && push!(sw.heap, (d, u, v))
    return nothing
end

# Crossings at the radius d where point id enters or leaves the segment at position t:
# two points entering (or leaving) at one end together, or a sphere appearing on
# another, meet then, but enter in their order after meeting, so no swap reports it.
# Walks from the neighbors l and r outward over the points at t (up to SWEEP_POS_SLACK)
# and records the monotonicity values of id with them that fall at d up to rounding (a
# superset is harmless: every value recorded is a critical value on the segment).
function sweep_coincident!(sw::MonoSweep, id::Int, d::Float64, (t, l, r), rec)
    for (w, step) in ((l, regress), (r, advance))
        while w != 0 && abs(pt_key(sw.clock.tab, w, d)[1] - t) <= SWEEP_POS_SLACK
            sweep_record_crossing!(sw, id, w, d, rec)
            w = sweep_neighbor(sw, w, step)
        end
    end
end

function sweep_record_crossing!(sw::MonoSweep, id::Int, w::Int, d::Float64, (out, lo, hi))
    k, l = minmax(pt_vertex(id), pt_vertex(w))
    k == l && return nothing
    tab = sw.clock.tab
    x = crossing(tab.L2, tab.t0[k], tab.h2[k], tab.t0[l], tab.h2[l])[2]
    (!isnan(x) && abs(x - d) <= SWEEP_COINCIDENT * d) && in_range_push!(out, x, lo, hi)
    return nothing
end

function sweep_insert!(sw::MonoSweep, id::Int, d::Float64, rec)
    sw.clock.delta = d
    # The exact entry position, not one recomputed with rounding; ties go by speed.
    t = pt_entry(sw.clock.tab, id)
    sweep_add!(sw, id, (t, pt_key(sw.clock.tab, id, d)[2]))
    l, r = sweep_neighbor(sw, id, regress), sweep_neighbor(sw, id, advance)
    sweep_coincident!(sw, id, d, (t, l, r), rec)
    l == 0 || sweep_schedule!(sw, l, id)
    r == 0 || sweep_schedule!(sw, id, r)
    return nothing
end

function sweep_remove!(sw::MonoSweep, id::Int, d::Float64, rec)
    sw.slot[id] == 0 && return nothing
    sw.clock.delta = d
    l, r = sweep_neighbor(sw, id, regress), sweep_neighbor(sw, id, advance)
    sweep_drop!(sw, id)
    # A point leaves at the end of the segment it moves toward.
    sweep_coincident!(sw, id, d, (pt_sign(id) < 0 ? 0.0 : 1.0, l, r), rec)
    (l != 0 && r != 0) && sweep_schedule!(sw, l, r)
    return nothing
end

# A scheduled swap; stale unless u and v are still adjacent in this order.
function sweep_swap!(sw::MonoSweep, (d, u, v), (out, lo, hi))
    (sw.slot[u] != 0 && sweep_neighbor(sw, u, advance) == v) || return nothing
    sw.clock.delta = max(sw.clock.delta, d)
    su, sv = sw.slot[u], sw.slot[v]
    sw.clock.item[su], sw.clock.item[sv] = v, u
    sw.slot[u], sw.slot[v] = sv, su
    tab = sw.clock.tab
    k, l = minmax(pt_vertex(u), pt_vertex(v))
    isnan(crossing(tab.L2, tab.t0[k], tab.h2[k], tab.t0[l], tab.h2[l])[2]) ||
        in_range_push!(out, d, lo, hi)
    l = sweep_neighbor(sw, v, regress)
    r = sweep_neighbor(sw, u, advance)
    l == 0 || sweep_schedule!(sw, l, v)
    r == 0 || sweep_schedule!(sw, u, r)
    return nothing
end

# Puts the points on the segment at radius lo into the set (sorted first, O(m log m)),
# and returns the later insertions and removals up to hi, sorted: (radius, 0 = insert /
# 2 = remove, point).
function sweep_init!(sw::MonoSweep, lo::Float64, hi::Float64)
    empty!(sw.set)
    fill!(sw.slot, 0)
    sw.free = collect(length(sw.slot):-1:1)
    sw.heap = BinaryMinHeap{Tuple{Float64,Int,Int}}()
    tab = sw.clock.tab
    sw.clock.delta = lo
    initial = Int[]
    events = Tuple{Float64,Int,Int}[]
    for id in eachindex(sw.slot)
        din, dout = pt_span(tab, id)
        (din <= dout && din <= hi && dout > lo) || continue
        din <= lo ? push!(initial, id) : push!(events, (din, 0, id))
        dout <= hi && push!(events, (dout, 2, id))
    end
    sort!(initial; by = id -> pt_key(tab, id, lo))
    for id in initial
        sweep_add!(sw, id, pt_key(tab, id, lo))
    end
    for j in 1:(length(initial)-1)
        sweep_schedule!(sw, initial[j], initial[j+1])
    end
    return sort!(events)
end

# Pushes the monotonicity values in [lo, hi] of the segment of sw.clock.tab onto out.
function mono_values!(out, sw::MonoSweep, lo::Float64, hi::Float64)
    # Events a little past hi still run: a removal there can record a coincident value
    # just below hi (sweep_coincident!).
    reach = hi + SWEEP_COINCIDENT * hi
    events = sweep_init!(sw, lo, reach)
    rec = (out, lo, hi)
    i = 1
    while i <= length(events) || !isempty(sw.heap)
        # At equal radii: insertions, then swaps, then removals.
        swap_first =
            !isempty(sw.heap) &&
            (i > length(events) || (first(sw.heap)[1], 1) < (events[i][1], events[i][2]))
        if swap_first
            e = pop!(sw.heap)
            e[1] > reach && break
            sweep_swap!(sw, e, rec)
        else
            d, kind, id = events[i]
            kind == 0 ? sweep_insert!(sw, id, d, rec) : sweep_remove!(sw, id, d, rec)
            i += 1
        end
    end
    return out
end

########################################################################
# Listing and sampling critical values.

# The vertex-vertex and vertex-segment values in [lo, hi] of the segment p[i]p[i+1]
# against the vertices of q (vertex-vertex for p[i], and for p[i+1] on the last one).
function edge_point_values!(out, tab::EdgeTable, last::Bool, lo, hi)
    for k in eachindex(tab.t0)
        in_range_push!(
            out,
            vertex_segment_value(tab.t0[k], tab.h2[k], tab.da2[k], tab.db2[k]),
            lo,
            hi,
        )
        in_range_push!(out, sqrt(tab.da2[k]), lo, hi)
        last && in_range_push!(out, sqrt(tab.db2[k]), lo, hi)
    end
end

# All monotonicity values in [lo, hi] of one segment, pair by pair: O(m^2).
function brute_mono_values!(out, tab::EdgeTable, lo, hi)
    m = length(tab.t0)
    for k in 1:m, l in (k+1):m
        _, d = crossing(tab.L2, tab.t0[k], tab.h2[k], tab.t0[l], tab.h2[l])
        isnan(d) || in_range_push!(out, d, lo, hi)
    end
end

# The critical values in [lo, hi] from the segments of one curve against the vertices
# of the other.
function list_critical!(out, s::SegTable, lo::Float64, hi::Float64, brute::Bool)
    tab = EdgeTable(other_count(s))
    sw = MonoSweep(tab)
    for i in eachindex(s.L2)
        edge_table!(tab, s, i)
        edge_point_values!(out, tab, i == length(s.L2), lo, hi)
        brute ? brute_mono_values!(out, tab, lo, hi) : mono_values!(out, sw, lo, hi)
    end
    return out
end

"""
    exact_critical_values(p, q, lo, hi; brute = false) -> Vector{Float64}
    exact_critical_values(t::ExactTables, lo, hi; brute = false)

The distinct critical values of the vertex vectors `p` and `q` (no consecutive
duplicates) in [lo, hi], sorted: by the kinetic sweep, or with `brute = true` by trying
every pair of vertices against every segment (O(n m (n + m))).
"""
exact_critical_values(p, q, lo::Float64, hi::Float64; brute::Bool = false) =
    exact_critical_values(exact_tables(p, q), lo, hi; brute)

function exact_critical_values(
    t::ExactTables,
    lo::Float64,
    hi::Float64;
    brute::Bool = false,
)
    out = Float64[]
    list_critical!(out, t.pq, lo, hi, brute)
    list_critical!(out, t.qp, lo, hi, brute)
    return unique!(sort!(out))
end

# Two distinct random indices k < l in 1:m.
function random_pair(rng, m::Int)
    k = rand(rng, 1:m)
    l = rand(rng, 1:(m-1))
    l += l >= k
    return minmax(k, l)
end

# A uniformly random monotonicity value of segment i of `s`, or NaN.
function sample_mono(rng, s::SegTable, i::Int)
    k, l = random_pair(rng, other_count(s))
    return crossing(
        s.L2[i],
        at(s, s.t0, k, i),
        at(s, s.h2, k, i),
        at(s, s.t0, l, i),
        at(s, s.h2, l, i),
    )[2]
end

# A uniformly random vertex-segment value of segment i of `s`.
function sample_vertex_segment(rng, s::SegTable, i::Int)
    k = rand(rng, 1:other_count(s))
    return vertex_segment_value(
        at(s, s.t0, k, i),
        at(s, s.h2, k, i),
        at(s, s.d2, k, i),
        at(s, s.d2, k, i + 1),
    )
end

# How many of `draws` critical values drawn uniformly from all of them fall on each
# kind: counts[1] vertex-vertex values; counts[2:5] per segment, vertex-segment values
# of the segments of p and of q, then monotonicity values of the segments of p and q.
# (Pairs of vertices with no monotonicity value count too; their draws come out NaN and
# are discarded, which keeps the rest uniform.)
function sample_counts(rng, t::ExactTables, draws::Int)
    n, m = space_size(t)
    w = (
        n * m,
        (n - 1) * m,
        (m - 1) * n,
        (n - 1) * binomial(m, 2),
        (m - 1) * binomial(n, 2),
    )
    total = sum(w)
    vv = 0
    counts = (zeros(Int, n - 1), zeros(Int, m - 1), zeros(Int, n - 1), zeros(Int, m - 1))
    for _ in 1:draws
        x = rand(rng, 1:total)
        kind = 1
        while x > w[kind]
            x -= w[kind]
            kind += 1
        end
        kind == 1 ? (vv += 1) : (c = counts[kind-1]; c[rand(rng, eachindex(c))] += 1)
    end
    return vv, counts
end

# 4 N^2 draws from the pairs of vertex and segment and the triples of two vertices and a
# segment that define the critical values, those in (lo, hi), sorted. A pair of vertices
# with no monotonicity value on the segment is drawn too and dropped: between two
# consecutive draws there are O(N log N) of these candidates with high probability, and
# so at most that many critical values, as the analysis needs. The draws are counted per segment first and then made segment by segment, so that each reads
# one column of the tables, in cache, rather than a random place in them: the same
# distribution as drawing them one by one.
function exact_sample(rng, t::ExactTables, lo::Float64, hi::Float64)
    n, m = space_size(t)
    vv, counts = sample_counts(rng, t, 4 * max(n, m)^2)
    R = Float64[]
    keep!(x) = (lo < x < hi) && push!(R, x)
    for _ in 1:vv
        keep!(sqrt(t.pq.d2[rand(rng, 1:m), rand(rng, 1:n)]))
    end
    sample_segments!(keep!, sample_vertex_segment, (rng, t.pq), counts[1])
    sample_segments!(keep!, sample_vertex_segment, (rng, t.qp), counts[2])
    sample_segments!(keep!, sample_mono, (rng, t.pq), counts[3])
    sample_segments!(keep!, sample_mono, (rng, t.qp), counts[4])
    return sort!(R)
end

# counts[i] draws of `draw` on segment i of s, for every segment.
function sample_segments!(keep!::K, draw::F, (rng, s), counts) where {K,F}
    for i in eachindex(counts), _ in 1:counts[i]
        keep!(draw(rng, s, i))
    end
end

########################################################################
# The search.

# The search works on `ctx = (tables, slack)`: the ExactTables of the two curves and the
# slack in ulps.
exact_accepts(ctx, c::Float64) = decide_free_space(ctx.tables, nextfloat(c, ctx.slack)^2)

# In the sorted vector R (with `lo` rejected and `hi` accepted around it): the largest
# rejected and smallest accepted values.
function exact_bracket(ctx, R::Vector{Float64}, lo::Float64, hi::Float64)
    i, j = 0, length(R) + 1
    while j - i > 1
        mid = (i + j) >>> 1
        exact_accepts(ctx, R[mid]) ? (j = mid) : (i = mid)
    end
    return (i == 0 ? lo : R[i]), (j > length(R) ? hi : R[j])
end

# The smallest critical value accepted by the decider, among those in (a, b] (b
# accepted, a rejected). The listing starts a little below a: the sweep starts from the
# order at its first radius, so a crossing there, computed an ulp above, can be missed.
function exact_search(ctx, a::Float64, b::Float64, brute::Bool)
    from = prevfloat(a, 4 * ctx.slack)
    S = filter!(x -> a < x < b, exact_critical_values(ctx.tables, from, b; brute))
    return exact_bracket(ctx, S, a, b)[2]
end


"""
    frechet_c_compute_exact(P, Q; rng = Random.default_rng(), method = :auto) -> Float64

The exact continuous Fréchet distance between the curves `P` and `Q`, by the
randomized algorithm of Har-Peled and Raichel ("The Fréchet distance revisited and
extended", Section 6): the distance is the smallest critical value of Alt and Godau
(vertex-vertex, vertex-segment and monotonicity events) that the exact decision
procedure (`frechet_c_decide_exact`) accepts. A random sample of 4 N² critical values
(N the larger number of vertices), sorted and binary searched with the decider, leaves
an interval that holds the distance and, with high probability, O(N log N) critical
values. These are listed (the monotonicity values by a kinetic sweep along each
segment, which keeps the points on the segment in a balanced search tree: O(log N) per
event), sorted, and searched again.

When the sample would cost more than about 35 decisions on the curves themselves
(N² ≥ 16 (n + m) d: long curves in low dimension), the interval comes instead from a
bisection with the fast decision procedure (`frechet_c_decide_exact` on `DecideCurve`s)
to a relative width of 1e-10; the critical values in it are computed from the
coordinates, segment by segment, and searched. With `method = :auto` (the default) this
is done without tables: the fast decision procedure also certifies the ends of the
interval, searches the critical values, and checks that nothing well below the result
is accepted; when a check fails, the tables are computed and the search is done on them.
With `method = :tables`, every decision of the search is made on the tables by the plain
cell by cell sweep, independent of the filters and the box search of the fast decision
procedure: the reference to check the default against. The two give the same result,
bit for bit, on all inputs tested.

**Slow:** without tables, O((n + m) d) memory and, besides the decisions, O(n m d) time
to compute the critical values near the distance (O(n m log N) when many points of a
segment are close together). With the tables: O(n m d) to tabulate, once, the distances
between vertices and segments (on all of Julia's threads; the result does not depend on
their number); then O(N² log N) to sample and sort, O(log N) decisions of O(n m) each,
and O(n m log N) to list (with high probability). The tables take 40 n m bytes and the
sample (drawn only when the bisection does not settle it) 32 N² bytes.
`frechet_c_compute` is usually much faster; this function is the reference it can be
checked against, with no tolerance.

The result is a critical value, computed in floating point, so exact up to its
rounding. Randomness affects only the running time: every `rng` gives the same result,
bit for bit. A candidate counts as accepted when the decider accepts it plus 64 ulps
(a relative 1.4e-14), which absorbs the different rounding of the decider and of the
critical value; two critical values closer than that cannot be told apart.

Takes any curves the other algorithms take (static or `Vector` points; matrix curves
are converted); curves of one vertex are allowed.
"""
function frechet_c_compute_exact(
    P::GPntSeq{N,T,PT},
    Q::GPntSeq{N,T,PT};
    rng::AbstractRNG = Random.default_rng(),
    method::Symbol = :auto,
)::Float64 where {N,T,PT<:AbstractVector{T}}
    who = "frechet_c_compute_exact"
    check_curves(P, Q, who; min_vertices = 1)
    method in (:auto, :tables) ||
        throw(ArgumentError("$who: method must be :auto or :tables (got :$method)"))
    P, Q = remove_consecutive_duplicates(P), remove_consecutive_duplicates(Q)
    if method == :auto && exact_bisects(P.pnts, Q.pnts)
        x = exact_without_tables(P, Q)
        x === nothing || return x
    end
    return exact_with_tables(P, Q, rng)
end

# The search on the tables of the two curves (no consecutive duplicates).
function exact_with_tables(P, Q, rng::AbstractRNG)
    p, q = P.pnts, Q.pnts
    # The tables take 40 n m bytes; the sample, drawn only without the bisection, 32 N^2
    # more (checked there).
    nmax = max(length(p), length(q))
    check_quadratic_size(length(p), length(q), 40, "frechet_c_compute_exact")
    tables = exact_tables(p, q)
    d2 = tables.pq.d2
    lo = sqrt(max(d2[1, 1], d2[end, end]))
    hi = sqrt(maximum(d2))
    (length(p) == 1 || length(q) == 1) && return hi
    ctx = (; tables, slack = EXACT_SLACK_ULPS)
    exact_accepts(ctx, lo) && return lo
    # The sample costs O(N^2) after the tables, a decision on the curves O((n + m) d):
    # bisect when the sample is the larger cost.
    if exact_bisects(p, q)
        x = exact_from_bisection(ctx, P, Q, lo, hi)
        x === nothing || return x
    end
    check_quadratic_size(nmax, nmax, 72, "frechet_c_compute_exact")
    a, b = exact_bracket(ctx, exact_sample(rng, tables, lo, hi), lo, hi)
    return exact_minimal(ctx, a, exact_search(ctx, a, b, false))
end

# Minimality: nothing well below x may be accepted, or the sweep missed a value; then
# the critical values in (a, x) are listed by brute force.
function exact_minimal(ctx, a::Float64, x::Float64)
    if decide_free_space(ctx.tables, prevfloat(x, ctx.slack)^2)
        Threads.atomic_add!(EXACT_FALLBACKS, 1)
        x = exact_search(ctx, a, x, true)
    end
    return x
end

# The relative width of the interval left by the bisection of exact_from_bisection.
const EXACT_BISECT_RTOL = 1.0e-10

# exact_from_bisection replaces the sample when N^2 >= EXACT_BISECT_FACTOR (n + m) d:
# about 35 decisions on the curves against 4 N^2 draws from the tables (measured: the
# bisection is faster on planar curves of 300 vertices and more, slower on curves of
# 100 to 200 vertices in R^963 and R^19200).
const EXACT_BISECT_FACTOR = 16

# The search without the sample: a bisection with the fast decision procedure
# (frechet_c_decide_exact on DecideCurves) narrows [lo, hi] to [a, b], b - a at most
# EXACT_BISECT_RTOL a. The decider of the search (on the tables) must reject a and accept
# b; then the critical values between a and `top`, b plus 16 times the slack in ulps
# (the accepted critical value can lie some ulps above b), are listed and searched as
# after the sample. The result is the smallest accepted critical value, as with the
# sample: the values up to a are rejected, the decider being monotone. When a or b is not
# certified, or no listed value is accepted, returns nothing, and the sample is drawn.
function exact_from_bisection(ctx, P, Q, lo::Float64, hi::Float64)
    a, b = exact_bisect(exact_decide_curves(P, Q)..., lo, hi)
    (exact_accepts(ctx, a) || !exact_accepts(ctx, b)) && return nothing
    S, top = narrow_candidates(P, Q, a, b, ctx.slack)
    x = exact_bracket(ctx, S, a, top)[2]
    x == top && return nothing
    return exact_minimal(ctx, a, x)
end
