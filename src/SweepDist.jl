###########################################################################
# SweepDist.jl
#
# Implementation of code to compute the sweeping distance. This is a
# similar entity to the continuous dynamic time warping distance, except
# that the price of matching two edges is approximated (fp_match_price:
# the average leash length times the length moved), and of course, we
# use our ve-r framework to compute/approximate it, together with the
# regular refinement for monotonicity, and refinement to get
# convergences.
###########################################################################

############################################################
# TreeVertex:
#    Encoding for a vertex in the VE Diagram how the search
#    procedure arrived ot it. I.e., encoding the reverse tree
#    that the algorithm computes.
#
# id: Vertex of the diagram
# val: The value associated with this vertex.
# id_prev: The vertex we arrvied from.
############################################################
@with_kw mutable struct TreeVertex
    id::Int64
    val::Float64 = 0.0
    #    r::Float64 = 0.0
    id_prev::Int64
end

function TreeVertex(_id::Int64)
    return TreeVertex(_id, 0.0, -1)
end





@with_kw mutable struct SDContext{N,T}
    P::PntSeq{N,T}
    Q::PntSeq{N,T}
    vertices::Dict{Int64,TreeVertex}
    handled::Dict{Int64,Bool}
    prev_map::Dict{Int64,Int64}
    pq::PriorityQueue{Int64,Float64,Base.Order.ForwardOrdering}
    n_p::Int64 = 0
    n_q::Int64 = 0
end

function SDContext(P::PntSeq{N,T}, Q::PntSeq{N,T}) where {N,T}
    return SDContext(
        P,
        Q,
        Dict{Int64,TreeVertex}(),
        Dict{Int64,Bool}(),
        Dict{Int64,Int64}(),
        PriorityQueue{Int64,Float64}(),
        cardin(P),
        cardin(Q),
    )
end


function get_event_points(c::SDContext{N,T}, id::Int64) where {N,T}
    P = c.P
    Q = c.Q
    i = EID_i(id)
    j = EID_j(id)
    if EID_i_is_vert(id)
        if EID_j_is_vert(id)
            return P[i], Q[j]
        end

        q = fp_nn_point(Q[j], Q[j+1], P[i])
        return P[i], q
    end

    if EID_j_is_vert(id)
        p = fp_nn_point(P[i], P[i+1], Q[j])
        return p, Q[j]
    end

    error("get_event_points: edge-edge events are not supported (event id $id)")
end



function SW_get_event_value(
    _id::Int64,
    _prev_id::Int64,
    prev_val::Float64,
    c::SDContext{N,T},
) where {N,T}
    o_p, o_q = get_event_points(c, _prev_id)
    n_p, n_q = get_event_points(c, _id)

    return prev_val + fp_match_price(o_p, n_p, o_q, n_q)
end


function SW_new_event(
    id::Int64,
    id_prev::Int64,
    value::Float64,
    c::SDContext{N,T},
) where {N,T}
    ev = TreeVertex(id)
    ev.id_prev = id_prev
    ev.val = value

    c.vertices[id] = ev

    return ev
end

function is_sw_schedule_event(
    dict::Dict{Int64,Bool},
    id::Int64,
    n_p::Int64,
    n_q::Int64,
)::Bool
    if haskey(dict, id)
        return false
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


function SW_schedule_event!(c::SDContext{N,T}, id::Int64, id_prev::Int64) where {N,T}
    if ! is_sw_schedule_event(c.handled, id, c.n_p, c.n_q)
        return
    end
    ev_prev::TreeVertex = c.vertices[id_prev]
    new_val = SW_get_event_value(id, id_prev, ev_prev.val, c)
    local ev
    if haskey(c.vertices, id)
        ev = c.vertices[id]
        if (ev.val < new_val)  # already better event scheduled.
            return
        end
        ev.val = new_val
    else
        ev = SW_new_event(id, id_prev, new_val, c)
    end

    ev.id_prev = id_prev
    if haskey(c.pq, id)
        delete!(c.pq, id)
    end
    enqueue!(c.pq, id, new_val)

    return ev
end



"""
    SweepDist_compute( P, Q )

A morphing between `P` and `Q` minimizing the approximate sweep price (the sum of
`fp_match_price` over matched edges), computed by a Dijkstra search over the
vertex-edge events. The minimum price is returned in `m.sol_value`; `m.leash` is
the longest leash of that morphing. The morphing need not be monotone.
"""
function SweepDist_compute(P::PntSeq{N,T}, Q::PntSeq{N,T}) where {N,T}
    check_curves(P, Q, "SweepDist_compute"; min_vertices = 2)
    c::SDContext{N,T} = SDContext(P, Q)

    start_id = EID(1, true, 1, true)

    id_end = EID(c.n_p, true, c.n_q, true)
    start_event = SW_new_event(start_id, start_id, 0.0, c)

    enqueue!(c.pq, start_id, start_event.val)

    end_event::TreeVertex = start_event
    iters = 0
    pq = c.pq
    while ! isempty(pq)
        tp = peek(pq)
        id = tp[1]
        value = tp[2]
        dequeue!(pq)

        if (haskey(c.handled, id))
            continue
        end

        ev::TreeVertex = c.vertices[id]
        c.handled[id] = true
        c.prev_map[id] = ev.id_prev
        iters = iters + 1


        i = EID_i(id)
        j = EID_j(id)

        if id == start_id
            SW_schedule_event!(c, EID(1, false, 1, true), id)
            SW_schedule_event!(c, EID(1, true, 1, false), id)
            continue
        end

        # Is it the *final* event?
        if (i == c.n_p) && (j == c.n_q)
            end_event = ev
            break
        end

        # Is it on the boundary of the final cell?
        if is_final_cell(ev.id, c.n_p, c.n_q)
            SW_schedule_event!(c, id_end, ev.id)
            continue
        end
        SW_schedule_event!(c, EID(i+1, true, j, false), id)
        SW_schedule_event!(c, EID(i, false, j+1, true), id)
    end

    pes, qes = f_r_extract_solution_basic(P, Q, id_end, c.prev_map)

    morph::Morphing{N,T} = Morphing_init(P, Q, pes, qes)
    morph.iters = iters
    morph.sol_value = end_event.val

    Morphing_verify_valid(morph)

    return morph
end


function SweepDist_compute_refine_mono(poly_a::PntSeq{N,T}, poly_b::PntSeq{N,T}) where {N,T}
    check_curves(poly_a, poly_b, "SweepDist_compute_refine_mono"; min_vertices = 2)
    ell = geom_length(poly_a) + geom_length(poly_b)

    P = poly_a
    Q = poly_b


    local m_new
    rounds::Int = 0
    while true
        rounds += 1
        check_rounds(rounds, "SweepDist_compute_refine_mono")
        m = SweepDist_compute(P, Q)
        if (Morphing_is_monotone(m))
            return m
        end

        m_new = Morphing_monotonize(m)
        rate = m_new.monotone_err / ell

        if (rate < 0.01)
            break
        end
        poly_a_2 = extract_refined_polygon(P, m.pes, 1)
        poly_b_2 = extract_refined_polygon(Q, m.qes, 1)

        P = poly_a_2
        Q = poly_b_2
    end

    return m_new
end

function SweepDist_compute_split(
    poly_a::PntSeq{N,T},
    poly_b::PntSeq{N,T},
    out::Vector{Morphing{N,T}},
    iters::Int64,
    max_vertices::Int64 = 1000000,
) where {N,T}
    check_curves(poly_a, poly_b, "SweepDist_compute_split"; min_vertices = 2)

    m = SweepDist_compute_refine_mono(poly_a, poly_b)
    Morphing_verify_valid(m)

    push!(out, deepcopy(m))
    P = m.P
    Q = m.Q
    for i in 1:iters
        P = split_edges(P)
        Q = split_edges(Q)
        mi = SweepDist_compute(P, Q)
        Morphing_verify_valid(mi)
        push!(out, deepcopy(mi))
        if ((cardin(P) + cardin(Q)) > max_vertices)
            break
        end
    end
end


###################################################################

"""
    SweepDist_lb_compute

Compute a lower bound on the Sweep distance between the two curves P
and Q. Current implementation is silly, and very similar to the
discrete Frechet calculation.

"""
function SweepDist_lb_compute(P::PntSeq{D,T}, Q::PntSeq{D,T}) where {D,T}
    check_curves(P, Q, "SweepDist_lb_compute"; min_vertices = 2)

    pq = PriorityQueue{Tuple{Int64,Int64},Float64}()

    n_p = cardin(P)
    n_q = cardin(Q)

    # Inf marks a cell not reached yet (0 is a legitimate value).
    dp = fill(Inf, n_p, n_q)
    dp_dec_i = falses(n_p, n_q)
    handled = falses(n_p, n_q)

    # edge_lb: for the grid edge (i, j) -> (i+1, j), the minimum of the distance
    # between the edge PA[i] PA[i+1] and the edges PB[j-1] PB[j] and PB[j] PB[j+1]
    # of the two adjacent cells: a lower bound on the leash along that edge.
    function edge_lb(PA::PntSeq{D,T}, i, PB::PntSeq{D,T}, j) where {D,T}
        lb::Float64 = typemax(Float64)
        if i == cardin(PA)
            return 0.0
        end
        if j > 1
            d = fp_dist_segment_segment(PA[i], PA[i+1], PB[j-1], PB[j])
            lb = min(lb, d)
        end
        if j < cardin(PB)
            d = fp_dist_segment_segment(PA[i], PA[i+1], PB[j], PB[j+1])
            lb = min(lb, d)
        end
        return lb
    end

    function edge_v(
        PA::PntSeq{D,T},
        i::Int64,
        PB::PntSeq{D,T},
        j::Int64,
    )::Float64 where {D,T}
        if (i == cardin(PA))
            return 0.0
        end
        lb = edge_lb(PA, i, PB, j)
        len = fp_dist(PA[i], PA[i+1])
        return len * lb
    end


    function push_value(i::Int64, j::Int64, ia::Int64, ja, new_val::Float64)
        #new_val = val + fp_dist( P[ ia ], Q[ ja ] )

        if dp[ia, ja] <= new_val
            return
        end
        dp[ia, ja] = new_val
        pq[ia, ja] = new_val  # Push to queue
        dp_dec_i[ia, ja] = (i < ia)
    end

    dp[1, 1] = 0.0
    enqueue!(pq, (1, 1), dp[1, 1])

    iters::Int64 = 0
    while ! isempty(pq)
        iters = iters + 1

        ele = peek(pq)
        i::Int64, j::Int64 = ele[1]
        value = ele[2]

        dequeue!(pq)
        if handled[i, j]
            continue
        end

        handled[i, j] = true
        if ((i == n_p) && (j == n_q))
            break
        end


        # Now we need to schedule the next two adjacent entries..
        if (i < n_p)
            push_value(i, j, i + 1, j, value + edge_v(P, i, Q, j))
        end
        if (j < n_q)
            push_value(i, j, i, j+1, value + edge_v(Q, j, P, i))
        end
    end

    m = d_frechet_extract_solution(P, Q, dp_dec_i, n_p, n_q)
    m.iters = iters

    m.sol_value = dp[n_p, n_q]

    return m
end

# SweepDist.jl end of file
#############################################################################
