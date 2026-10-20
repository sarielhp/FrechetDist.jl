# Originally contributed by S. Har-Peled
# under MIT License

#-------------------------------------------------------------
# Frechet
#
# A library to compute various versions of the Frechet distance.
#-------------------------------------------------------------

const DictVERType = Dict{Int64,Int64};
#DictHandledType = Dict{Int64, Bool};


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

#DictVERType = Dict{Int64, TreeVertex};
#HeapVERType = BinaryMinHeap{TreeVertex};

const HeapVerticesType = BinaryMinHeap{Tuple{Float64,Int64}};


#@with_kw
mutable struct FRContext{N,T}
    P::PntSeq{N,T}
    Q::PntSeq{N,T}
    p_offs::Vector{Float64}
    q_offs::Vector{Float64}
    #    handled::DictHandledType;
    dict::DictVERType
    heapAlt::HeapVerticesType
    f_offsets::Bool
    n_p::Int64 #0
    n_q::Int64 #0
    f_upper_bound::Bool
    upper_bound::T
end

function FR_Context(P::PntSeq{N,T}, Q::PntSeq{N,T}) where {N,T}
    #    d_h = DictHandledType();
    d_ver = DictVERType()
    h = HeapVerticesType()
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
            #                println( "Larger than 0!" );
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
    #=
    if  ( any( isnan, p ) )
        println( seg );
        println( qr );
        @assert( false );
    end
    =#

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
        t = 0
        p = R[i]
    end
    if (t < 1.0) && ((1.0 - t) * len <= noise)
        t = 1.0
        p = R[i+1]
    end
    #        println( "TTT= ", t );

    @assert(! any(isnan, p))
    return EventPoint(p, Int64(i), PT_ON_EDGE, t, 0.0)
end


@inline function is_schedule_event(
    dict::DictVERType,
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
    dict::DictVERType,
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
    dict::DictVERType,
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

    if (length(rng) <= 2)
        return 0
    end

    s_p = P[first(rng)]
    s_q = P[last(rng)]

    seg_len_sq = fp_dist_sq(s_p, s_q)
    s_vec = s_q - s_p

    t::Float64 = 0
    curr::Point{N,T} = s_p
    leash::Float64 = 0
    for i = (first(rng)+1):(last(rng)-1)
        q, new_t = fp_nn_point_ext_ext(s_p, s_q, P[i], seg_len_sq, s_vec)
        #q = fp_nn_point( seg, P[ i ] );
        #new_t = fp_convex_coef( seg, q );
        if (new_t > t)
            t = new_t
            curr = q
        end
        leash = max(leash, fp_dist_sq(curr, P[i]))
    end
    return sqrt(leash)
end

eachindexButLast(x) = firstindex(x):(lastindex(x)-1)

#####################################################################
"""
    frechet_offests

Given P, and subset of indices of P defining a subpolygon, computes a
-approximation to the Frechet distance between every edge of the
subpolygon they define, and corresponding portion of P.  A fast, and
probably decent approximation to the optimal morhping. Finally, set
for every vertex of the subcurve, the frechet distance of this
subcurve. Thus, providing a fast estimates of the offsets.

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
