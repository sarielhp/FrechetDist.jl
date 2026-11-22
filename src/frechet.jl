# Originally contributed by S. Har-Peled
# under MIT License

#-------------------------------------------------------------
# Frechet
#
# A library to compute various versions of the Frechet distance.
#-------------------------------------------------------------



"""
    EID

A function that encodes a matching "edge" betwee two polygons.
"""
function EID(
    i::Int64 = 0,
    i_is_vert::Bool = false,
    j::Int64 = 0,
    j_is_vert::Bool = false,
)::Int64
    #---
    val::Int64 = (
        (i << 32) | (j << 2) | (convert(Int64, i_is_vert) << 1) |
        (convert(Int64, j_is_vert))
    )
    return val
end

@inline function EID_i(id::Int64)::Int64
    return (id >> 32)
end

@inline function EID_j(id::Int64)::Int64
    return (id & 0xfffffffc) >> 2
end

@inline function EID_ij_status(id::Int64)::Int64
    return id & 0x3
end

@inline function EID_i_is_vert(id::Int64)::Bool
    return ((id & 0x2) != 0)
end

@inline function EID_j_is_vert(id::Int64)::Bool
    return ((id & 0x1) != 0)
end




"""
    PrevMap

The predecessor map of the VE search: event id => id of the event it was reached
from. It starts as a `Dict`, which is cheap when the search explores few events
(similar curves), and moves once to a dense array indexed by (i, j, status) when it
holds more than `limit` entries, which is faster and smaller when the search explores
most of the 4 n_p n_q events. Only lookups are used, so the representation does not
affect results.
"""
mutable struct PrevMap
    dict::Dict{Int64,Int64}
    dense::Vector{Int64}   # empty until the move; 0 marks a missing entry
    n_p::Int64
    n_q::Int64
    limit::Int64
end

"Smallest number of entries before a `PrevMap` turns dense. A `Ref` so tests can lower it."
const PREVMAP_MIN_DENSE = Ref(1024)

function PrevMap(n_p::Int64, n_q::Int64)
    # Stay a Dict when the dense array would need more than a quarter of the memory.
    dense_ok = 32.0 * n_p * n_q <= Sys.total_memory() / 4
    limit = dense_ok ? max(PREVMAP_MIN_DENSE[], (n_p * n_q) >> 4) : typemax(Int64)
    return PrevMap(Dict{Int64,Int64}(), Int64[], n_p, n_q, limit)
end

@inline pm_slot(m::PrevMap, id::Int64) =
    ((EID_j(id) - 1) * m.n_p + (EID_i(id) - 1)) * 4 + EID_ij_status(id) + 1
@inline pm_in_range(m::PrevMap, id::Int64) =
    (1 <= EID_i(id) <= m.n_p) && (1 <= EID_j(id) <= m.n_q)

@inline function Base.haskey(m::PrevMap, id::Int64)
    isempty(m.dense) && return haskey(m.dict, id)
    return pm_in_range(m, id) && (m.dense[pm_slot(m, id)] != 0)
end

@inline function Base.getindex(m::PrevMap, id::Int64)
    isempty(m.dense) && return m.dict[id]
    v = pm_in_range(m, id) ? m.dense[pm_slot(m, id)] : 0
    v == 0 && throw(KeyError(id))
    return v
end

function Base.setindex!(m::PrevMap, v::Int64, id::Int64)
    if !isempty(m.dense)
        m.dense[pm_slot(m, id)] = v
        return m
    end
    m.dict[id] = v
    length(m.dict) > m.limit && pm_make_dense!(m)
    return m
end

function pm_make_dense!(m::PrevMap)
    m.dense = zeros(Int64, 4 * m.n_p * m.n_q)
    for (k, v) in m.dict
        m.dense[pm_slot(m, k)] = v
    end
    m.dict = Dict{Int64,Int64}()
    return m
end

"A predecessor map: the VE search uses a `PrevMap`, the sweep search a `Dict`."
const PrevLookup = Union{Dict{Int64,Int64},PrevMap}

mutable struct FRContext{N,T}
    P::PntSeq{N,T}
    Q::PntSeq{N,T}
    p_offs::Vector{Float64}
    q_offs::Vector{Float64}
    dict::PrevMap
    heapAlt::BinaryMinHeap{Tuple{Float64,Int64}}
    f_offsets::Bool
    n_p::Int64 #0
    n_q::Int64 #0
    f_upper_bound::Bool
    upper_bound::T
end

function FR_Context(P::PntSeq{N,T}, Q::PntSeq{N,T}) where {N,T}
    d_ver = PrevMap(Int64(cardin(P)), Int64(cardin(Q)))
    h = BinaryMinHeap{Tuple{Float64,Int64}}()
    return FRContext(
        P,
        Q,
        Vector{Float64}(),
        Vector{Float64}(),
        d_ver,
        h,
        false,
        Int64(cardin(P)),
        Int64(cardin(Q)),
        false,
        zero(T),  # Upper bound
    )
end


function ve_event_value(c::FRContext{N,T}, id::Int64) where {N,T}
    P = c.P
    Q = c.Q
    i = EID_i(id)
    j = EID_j(id)
    penalty::T = 0

    if EID_i_is_vert(id)
        if EID_j_is_vert(id)
            d = fp_dist(P[i], Q[j])
            # Old code: Bug? Missing if check
            if c.f_offsets
                penalty = c.p_offs[i] + c.q_offs[j]
                return d - penalty, penalty
            end
            return fp_dist(P[i], Q[j]), 0.0
        end
        d = fp_dist_point_segment(c.P[i], Q[j], Q[j+1])
        if c.f_offsets
            #            if  ( c.q_offs[ j ] > 0 )
            #            end
            #--------------------------------------------------------------
            # Old code: Bug?
            #            return  d - c.q_offs[ j ];
            penalty = max(c.q_offs[j], c.q_offs[j+1]) + c.p_offs[i]
            return d - penalty, penalty
        end
        return d, 0.0
    end

    if EID_j_is_vert(id)
        d = fp_dist_point_segment(Q[j], P[i], P[i+1])
        if c.f_offsets
            penalty = max(c.p_offs[i], c.p_offs[i+1]) + c.q_offs[j]
            return d - penalty, penalty
        end
        return d, 0.0
    end

    error("ve_event_value: edge-edge events are not supported (event id $id)")
end


@inline function is_start_event(id::Int64)
    return ((EID_i(id) == 1) && (EID_j(id) == 1) && EID_i_is_vert(id) && EID_j_is_vert(id))
end

@inline function is_final_cell(id::Int64, n_p::Int64, n_q::Int64)
    return ((EID_i(id) == (n_p - 1)) && (EID_j(id) == (n_q - 1)))
end

function f_r_create_event(
    R::PntSeq{N,T},
    i::Integer,
    is_vert::Bool,
    qr::Point{N,T},
) where {N,T}
    if (is_vert)
        ev = EventPoint(R[i], Int64(i), PT_VERTEX, 0.0, 0.0)
        return ev
    end

    p::Point{N,T}, t::T = fp_nn_point_ext(R[i], R[i+1], qr)


    #seg = Segment( R[ i ], R[ i + 1 ] );
    #p = fp_nn_point( seg, qr );

    #t::Float64 = fp_convex_coef( seg, p );

    (0.0 <= t <= 1.0) || error(
        "f_r_create_event: the nearest-point parameter t = $t on edge $i is outside " *
        "[0, 1] (query point $qr)",
    )
    # Snap to the vertices only when the displacement is rounding noise for the
    # coordinates involved (independent of the edge length).
    noise = 4 * eps(Float64) * max(maximum(abs, R[i]), maximum(abs, R[i+1]))
    len = fp_dist(R[i], R[i+1])
    if (0 < t) && (t * len <= noise)
        t = 0.0
        p = R[i]
    end
    if (t < 1.0) && ((1.0 - t) * len <= noise)
        t = 1.0
        p = R[i+1]
    end

    @assert(! any(isnan, p))
    return EventPoint(p, Int64(i), PT_ON_EDGE, t, 0.0)
end


@inline function is_schedule_event(
    dict::PrevLookup,
    id::Int64,
    n_p::Int64,
    n_q::Int64,
)::Bool
    if haskey(dict, id)
        return false
    end
    if (EID_i(id) < n_p) && (EID_j(id) < n_q)
        return true
    end
    if (EID_i(id) > n_p) || (EID_j(id) > n_q)
        return false
    end
    if (EID_i(id) >= n_p) && (! EID_i_is_vert(id))
        return false
    end
    if (EID_j(id) >= n_q) && (! EID_j_is_vert(id))
        return false
    end
    return true
end

@inline function f_r_schedule_event(
    id::Int64,
    prev_id::Int64,
    c::FRContext{N,T},
) where {N,T}
    if ! is_schedule_event(c.dict, id, c.n_p, c.n_q)
        return
    end
    new_val, penalty = ve_event_value(c, id)
    if (c.f_upper_bound && (c.upper_bound < new_val))
        return
    end
    ### Also... Blocks the event from being considered again...
    c.dict[id] = prev_id
    push!(c.heapAlt, (new_val, id))
end


function f_r_extract_solution_basic(
    P::PntSeq{N,T},
    Q::PntSeq{N,T},
    end_event_id::Int64,
    dict::PrevLookup,
) where {N,T}
    #############################################################3
    # Extracting the solution
    #
    # The final event is always the same - final vertex of P versus
    # final vertex of Q.
    pes = Vector{EventPoint{N,T}}()
    qes = Vector{EventPoint{N,T}}()

    idx = end_event_id
    pex::EventPoint = f_r_create_event(P, EID_i(idx), EID_i_is_vert(idx), Q[EID_j(idx)])
    qex::EventPoint = f_r_create_event(Q, EID_j(idx), EID_j_is_vert(idx), P[EID_i(idx)])

    push!(pes, pex)
    push!(qes, qex)

    curr = end_event_id
    while ! is_start_event(curr)
        prev = curr
        curr = dict[prev]

        #        id = curr.id;
        id_i = EID_i(curr)
        id_j = EID_j(curr)

        pe::EventPoint = f_r_create_event(P, id_i, EID_i_is_vert(curr), Q[id_j])
        qe::EventPoint = f_r_create_event(Q, id_j, EID_j_is_vert(curr), P[id_i])
        push!(pes, pe)
        push!(qes, qe)
    end

    reverse!(pes)
    reverse!(qes)

    return pes, qes
end

function f_r_extract_solution(
    P::PntSeq{N,T},
    Q::PntSeq{N,T},
    end_event_id::Int64,
    dict::PrevLookup,
    c::FRContext{N,T},
) where {N,T}
    #############################################################3
    # Extracting the solution
    #
    # The final event is always the same - final vertex of P versus
    # final vertex of Q.
    pes = Vector{EventPoint{N,T}}()
    qes = Vector{EventPoint{N,T}}()

    idx = end_event_id
    pex::EventPoint = f_r_create_event(P, EID_i(idx), EID_i_is_vert(idx), Q[EID_j(idx)])
    qex::EventPoint = f_r_create_event(Q, EID_j(idx), EID_j_is_vert(idx), P[EID_i(idx)])

    if (c.f_offsets)
        val, penalty = ve_event_value(c, idx)
        pex = ev_update_penalty(pex, penalty)
        #pex.penalty = penalty;
    end
    push!(pes, pex)
    push!(qes, qex)

    curr = end_event_id
    while ! is_start_event(curr)
        prev = curr
        curr = dict[prev]

        #        id = curr.id;
        id_i = EID_i(curr)
        id_j = EID_j(curr)

        pe::EventPoint = f_r_create_event(P, id_i, EID_i_is_vert(curr), Q[id_j])
        qe::EventPoint = f_r_create_event(Q, id_j, EID_j_is_vert(curr), P[id_i])
        if (c.f_offsets)
            val, penalty = ve_event_value(c, curr)
            pe = ev_update_penalty(pe, penalty)
        end
        push!(pes, pe)
        push!(qes, qe)
    end

    reverse!(pes)
    reverse!(qes)

    return pes, qes
end



function eid_same_status(arr::Vector{Int64}, curr::Int64, len::Int64)
    eid = arr[curr]
    status = EID_ij_status(arr[curr])
    t = curr + 1
    while (t <= len)
        if EID_ij_status(arr[t]) != status
            break
        end
        t = t + 1
    end
    return t - 1
end


function f_r_extract_solution_ids(
    P::PntSeq{N,T},
    Q::PntSeq{N,T},
    end_event_id::Int64,
    start_event_id::Int64,
    dict,
) where {N,T}
    out_arr = Vector{Int64}()

    push!(out_arr, end_event_id)

    curr = end_event_id
    while (curr != start_event_id)
        push!(out_arr, curr)
        prev = curr
        haskey(dict, prev) || error(
            "f_r_extract_solution_ids: event $prev has no predecessor; the search did " *
            "not reach the end event (an upper bound below the VE distance?)",
        )
        curr = dict[prev]
    end
    push!(out_arr, curr)

    reverse!(out_arr)

    return out_arr
end

"""
    max_leash

    Gets a segment on one polygon s_a s_b, and a chain P[low:hi]. Updates the
    min/max estimates from the leash length. Lower bound is simply the distance
    to nearest point on th esegment. The max is the result of brute force
    monotonization.
"""
@inline function max_leash(
    bounds::Tuple{T,T},
    s_a::Point{N,T},
    s_b::Point{N,T},
    P::PntSeq{N,T},
    idx::UnitRange{Int64},
) where {N,T}
    l_min, l_max = bounds
    #seg = Segment( s_a, s_b );
    len_seg_sq::Float64 = fp_dist_sq(s_a, s_b)
    if (len_seg_sq == 0.0)
        for j in idx
            l_min = max(l_min, fp_dist(s_a, P[j]))
        end
        return l_min, max(l_min, l_max)
    end

    max_t = 0.0
    max_s = s_a

    seg_len_sq = fp_dist_sq(s_a, s_b)  # Squared distance between s_a and s_b
    s_vec = s_b - s_a

    for j in idx
        p = P[j]

        sq, new_t = fp_nn_point_ext_ext(s_a, s_b, p, seg_len_sq, s_vec)

        dst = fp_dist(sq, p)
        l_min = max(l_min, dst)
        if (new_t >= max_t)
            max_t = new_t
            max_s = sq #fp_convex_comb( p_a, p_b, max_t );
        end

        l_max = max(l_max, fp_dist(max_s, p))
    end

    return l_min, l_max
end

# Processes the run of consecutive ids starting at arr[curr] that are
# vertex-on-segment events. V is the curve contributing the vertices, S the
# curve contributing the segment; vert_is_i tells whether V is the first curve.
# Returns the updated (l_min, l_max) and the next position in arr.
function leash_vertex_run(
    bounds::Tuple{T,T},
    arr::Vector{Int64},
    curr::Int64,
    V::PntSeq{N,T},
    S::PntSeq{N,T},
) where {N,T}
    l_min, l_max = bounds
    vert_is_i = EID_i_is_vert(arr[curr])
    low = curr
    hi = eid_same_status(arr, curr, length(arr))
    if low > hi
        return l_min, l_max, curr + 1
    end
    eid_start = arr[low]
    eid_end = arr[hi]
    v_start = vert_is_i ? EID_i(eid_start) : EID_j(eid_start)
    s_idx = vert_is_i ? EID_j(eid_start) : EID_i(eid_start)
    s_a = S[s_idx]
    s_b = S[s_idx+1]

    if low == hi
        return max(l_min, fp_dist_point_segment(V[v_start], s_a, s_b)), l_max, curr + 1
    end

    v_end = vert_is_i ? EID_i(eid_end) : EID_j(eid_end)
    l_min, l_max = max_leash((l_min, l_max), s_a, s_b, V, v_start:v_end)
    return l_min, l_max, hi + 1
end


function compute_leash_from_arr(
    P::PntSeq{N,T},
    Q::PntSeq{N,T},
    arr::Vector{Int64},
) where {N,T}
    curr::Int64 = 1
    len::Int64 = length(arr)
    l_max::T = 0.0
    l_min::T = 0.0

    while curr <= len
        eid = arr[curr]
        i_vert = EID_i_is_vert(eid)
        j_vert = EID_j_is_vert(eid)

        if i_vert && j_vert
            l_min = max(l_min, fp_dist(P[EID_i(eid)], Q[EID_j(eid)]))
            curr = curr + 1
        elseif i_vert
            l_min, l_max, curr = leash_vertex_run((l_min, l_max), arr, curr, P, Q)
        elseif j_vert
            l_min, l_max, curr = leash_vertex_run((l_min, l_max), arr, curr, Q, P)
        else
            @assert(false)
        end
    end

    l_max = max(l_min, l_max)
    return l_min, l_max
end


"""
    frechet_ve_r_compute_range

    Return a range of distances. The lower end is the VE-Frechet
    distance, the upper end is a constructive upper bound on the
    Frechet distance.

"""
function frechet_ve_r_compute_range(
    P::PntSeq{N,T},
    Q::PntSeq{N,T},
    upper_bound::T,
) where {N,T}
    check_curves(P, Q, "frechet_ve_r_compute_range"; min_vertices = 2)
    check_nonnegative(upper_bound, "frechet_ve_r_compute_range", "upper_bound")
    c::FRContext{N,T} = FR_Context(P, Q)
    c.f_upper_bound = true
    c.upper_bound = upper_bound

    n_pm::Int64 = c.n_p - 1
    n_qm::Int64 = c.n_q - 1

    start_id = EID(1, true, 1, true)

    end_id = EID(c.n_p, true, c.n_q, true)
    f_r_schedule_event(start_id, start_id, c)

    heapAlt = c.heapAlt
    while ! isempty(heapAlt)
        tp = pop!(heapAlt)
        id = last(tp)
        value = first(tp)

        i = EID_i(id)
        j = EID_j(id)

        if id == start_id
            f_r_schedule_event(EID(1, false, 1, true), id, c)
            f_r_schedule_event(EID(1, true, 1, false), id, c)
            continue
        end

        if (i >= n_pm) && (j >= n_qm)
            # Is it the *final* event?
            if (i == c.n_p) && (j == c.n_q)
                break
            end
            if is_final_cell(id, c.n_p, c.n_q)
                f_r_schedule_event(end_id, id, c)
                continue
            end
        end
        f_r_schedule_event(EID(i+1, true, j, false), id, c)
        f_r_schedule_event(EID(i, false, j+1, true), id, c)
    end

    if ! haskey(c.dict, end_id)
        # The bound is below the VE distance (or equal to it up to rounding), so the
        # search was pruned before reaching the end. Compute the range without a bound,
        # as FEVER_compute_range does.
        m = frechet_ve_r_compute(P, Q)
        return m.leash, Morphing_monotonize(m).leash
    end

    out_arr = f_r_extract_solution_ids(P, Q, end_id, start_id, c.dict)

    l_min, l_max = compute_leash_from_arr(P, Q, out_arr)

    return l_min, l_max
end



##########################################################################
##########################################################################
# Compute retractable Frechet. Essentially Prim/Dijkstra algorithm
#
# Returns a Morphing that encodes the solution
##########################################################################
function frechet_ve_r_compute_ext(
    P::PntSeq{N,T},
    Q::PntSeq{N,T},
    p_offs::Vector{Float64},
    q_offs::Vector{Float64},
    f_use_offsets::Bool = false,
) where {N,T}
    check_curves(P, Q, "frechet_ve_r_compute_ext"; min_vertices = 2)
    if f_use_offsets && (length(p_offs) != cardin(P) || length(q_offs) != cardin(Q))
        throw(
            ArgumentError(
                "frechet_ve_r_compute_ext: p_offs and q_offs need one entry per vertex of their curve",
            ),
        )
    end
    if f_use_offsets
        check_offsets(p_offs, "frechet_ve_r_compute_ext", "p_offs")
        check_offsets(q_offs, "frechet_ve_r_compute_ext", "q_offs")
    end
    c = FR_Context(P, Q)

    if f_use_offsets
        c.p_offs = p_offs
        c.q_offs = q_offs
        c.f_offsets = true
    end

    start_id = EID(1, true, 1, true)

    end_id = EID(c.n_p, true, c.n_q, true)
    f_r_schedule_event(start_id, start_id, c)

    #end_event_i::TreeVertex = start_event;
    iters = 0
    heapAlt = c.heapAlt
    while ! isempty(heapAlt)
        tp = pop!(heapAlt)
        id = last(tp)
        value = first(tp)
        iters = iters + 1

        i = EID_i(id)
        j = EID_j(id)

        if id == start_id #is_start_event( id )
            f_r_schedule_event(EID(1, false, 1, true), id, c)
            f_r_schedule_event(EID(1, true, 1, false), id, c)
            continue
        end

        # Is it the *final* event?
        if (i == c.n_p) && (j == c.n_q)
            break
        end

        # Is it on the boundary of the final cell?
        if is_final_cell(id, c.n_p, c.n_q)
            f_r_schedule_event(end_id, id, c)
            continue
        end
        f_r_schedule_event(EID(i+1, true, j, false), id, c)
        f_r_schedule_event(EID(i, false, j+1, true), id, c)
    end

    #@time
    #XXX
    pes, qes = f_r_extract_solution(P, Q, end_id, c.dict, c)

    morph::Morphing{N,T} = Morphing_init(P, Q, pes, qes)
    morph.iters = iters

    return morph
end


function frechet_ve_r_compute(P::PntSeq{N,T}, Q::PntSeq{N,T}) where {N,T}
    check_curves(P, Q, "frechet_ve_r_compute"; min_vertices = 2)
    return frechet_ve_r_compute_ext(P, Q, Vector{Float64}(), Vector{Float64}(), false)
end




"""
    frechet_width_approx

 2-approximation to the Frechet distance between
   seg = P[first(rng)]-P[last(rng)] and he polygon
   P[rng]
 Here, rng is a range i:j

This function implements a greedy matching alnog the segment - you
move to the cloest point on seg to the current point, making sure
never moving back.

"""
function frechet_width_approx(P::PntSeq{N,T}, rng::UnitRange{Int64} = 0:0)::T where {N,T}
    card = cardin(P)
    if (card <= 2)
        return 0
    end

    if (rng == 0:0)
        rng = 1:card
    end
    (isempty(rng) || (1 <= first(rng) && last(rng) <= card)) || throw(
        ArgumentError("frechet_width_approx: range $rng is outside the vertices 1:$card"),
    )

    if (length(rng) <= 2)
        return 0
    end

    return sqrt(width_leash_sq(P, rng, Inf))
end

# The squared greedy width of P[rng] against its spine (rng has at least 3 vertices).
# Stops early once sqrt of the running maximum exceeds `w`: the maximum only grows and
# sqrt is monotone, so the comparison with `w` is decided at that point.
function width_leash_sq(P::PntSeq{N,T}, rng::UnitRange{Int64}, w::Float64) where {N,T}
    s_p = P[first(rng)]
    s_q = P[last(rng)]

    seg_len_sq = fp_dist_sq(s_p, s_q)
    s_vec = s_q - s_p

    t::Float64 = 0
    curr::Point{N,T} = s_p
    leash::Float64 = 0
    for i = (first(rng)+1):(last(rng)-1)
        q, new_t = fp_nn_point_ext_ext(s_p, s_q, P[i], seg_len_sq, s_vec)
        if (new_t > t)
            t = new_t
            curr = q
        end
        leash = max(leash, fp_dist_sq(curr, P[i]))
        (sqrt(leash) > w) && return leash
    end
    return leash
end

"""
    frechet_width_exceeds( P, rng, w )

Whether `frechet_width_approx(P, rng) > w`, computed with the same arithmetic but
stopping as soon as the answer is known.
"""
function frechet_width_exceeds(
    P::PntSeq{N,T},
    rng::UnitRange{Int64},
    w::Float64,
) where {N,T}
    (cardin(P) <= 2 || length(rng) <= 2) && return zero(T) > w
    return sqrt(width_leash_sq(P, rng, w)) > w
end

eachindexButLast(x) = firstindex(x):(lastindex(x)-1)

#####################################################################
"""
    frechet_offsets( P, p_indices )

Given `P` and the increasing indices `p_indices` of a subcurve of `P`, returns, for
every edge `p_indices[i] p_indices[i+1]` of the subcurve, a fast upper estimate
(`frechet_width_approx`) of the Frechet distance between that edge and the portion
of `P` it replaces: one value per edge.

"""
function frechet_offsets(P::PntSeq{N,T}, p_indices::Vector{Int64}) where {N,T}
    check_indices(p_indices, cardin(P), "frechet_offsets")
    offs = Vector{Float64}()

    for i in eachindexButLast(p_indices)
        d = frechet_width_approx(P, p_indices[i]:p_indices[i+1])
        push!(offs, d)
    end
    return offs
end
