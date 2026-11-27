###########################################################################
# SweepDist.jl
#
# Code to compute the sweeping distance. This is similar to the
# continuous dynamic time warping distance, except
# that the price of matching two edges is approximated (fp_match_price:
# the average leash length times the length moved), and of course, we
# use our ve-r framework to compute/approximate it, together with the
# regular refinement for monotonicity, and refinement to get
# convergence.
###########################################################################

# The events of the sweep search, by event id: the best value found so far and the
# event it was reached from (both valid once the event is created), and whether it
# has been popped (handled). Like PrevMap, it starts with Dicts, which are cheap when
# the search reaches few events (similar curves), and moves once to dense arrays
# indexed like PrevMap's when it holds more than `limit` events; the search on
# unrelated curves reaches almost all of them, and the arrays avoid the hashing.
mutable struct SweepStore
    dense::Bool
    n_p::Int64
    n_q::Int64
    limit::Int64
    vals::Vector{Float64}  # dense: allocated uninitialized; read only if created
    prev::Vector{Int64}
    created::BitVector
    handled::BitVector
    dvals::Dict{Int64,Tuple{Float64,Int64}}  # sparse: (value, prev)
    dhandled::Set{Int64}
end

function SweepStore(n_p::Int64, n_q::Int64)
    # Stay sparse when the arrays would need more than a quarter of the memory.
    dense_ok = 64.0 * n_p * n_q <= Sys.total_memory() / 4
    limit = dense_ok ? max(PREVMAP_MIN_DENSE[], (n_p * n_q) >> 4) : typemax(Int64)
    return SweepStore(
        false,
        n_p,
        n_q,
        limit,
        Float64[],
        Int64[],
        falses(0),
        falses(0),
        Dict{Int64,Tuple{Float64,Int64}}(),
        Set{Int64}(),
    )
end

function sw_make_dense!(st::SweepStore)
    slots = 4 * st.n_p * st.n_q
    st.vals = Vector{Float64}(undef, slots)
    st.prev = Vector{Int64}(undef, slots)
    st.created = falses(slots)
    st.handled = falses(slots)
    st.dense = true
    for (id, (v, p)) in st.dvals
        k = sw_slot(st, id)
        st.vals[k] = v
        st.prev[k] = p
        st.created[k] = true
    end
    for id in st.dhandled
        st.handled[sw_slot(st, id)] = true
    end
    st.dvals = Dict{Int64,Tuple{Float64,Int64}}()
    st.dhandled = Set{Int64}()
    return nothing
end

# The slot of an event with 1 <= i <= n_p and 1 <= j <= n_q (as pm_slot).
@inline sw_slot(st::SweepStore, id::Int64) =
    ((EID_j(id) - 1) * st.n_p + (EID_i(id) - 1)) * 4 + EID_ij_status(id) + 1

# The value of a created event, or Inf (larger than every price) if there is none.
@inline function sw_val_or_inf(st::SweepStore, id::Int64)
    st.dense || return first(get(st.dvals, id, (Inf, 0)))
    k = sw_slot(st, id)
    return st.created[k] ? st.vals[k] : Inf
end
@inline sw_val(st::SweepStore, id::Int64) =
    st.dense ? st.vals[sw_slot(st, id)] : first(st.dvals[id])
@inline sw_prev(st::SweepStore, id::Int64) =
    st.dense ? st.prev[sw_slot(st, id)] : last(st.dvals[id])
@inline sw_handled(st::SweepStore, id::Int64) =
    st.dense ? st.handled[sw_slot(st, id)] : (id in st.dhandled)

@inline function sw_set!(st::SweepStore, id::Int64, val::Float64, prev::Int64)
    if st.dense
        k = sw_slot(st, id)
        st.vals[k] = val
        st.prev[k] = prev
        st.created[k] = true
    else
        st.dvals[id] = (val, prev)
        length(st.dvals) > st.limit && sw_make_dense!(st)
    end
    return nothing
end

@inline function sw_handle!(st::SweepStore, id::Int64)
    st.dense ? (st.handled[sw_slot(st, id)] = true) : push!(st.dhandled, id)
    return nothing
end


struct SDContext{N,T,PT<:AbstractVector{T}}
    P::GPntSeq{N,T,PT}
    Q::GPntSeq{N,T,PT}
    store::SweepStore
    prev_map::PrevMap
    pq::PriorityQueue{Int64,Float64,Base.Order.ForwardOrdering}
    n_p::Int64
    n_q::Int64
end

function SDContext(P::GPntSeq{N,T,PT}, Q::GPntSeq{N,T,PT}) where {N,T,PT<:AbstractVector{T}}
    n_p, n_q = cardin(P), cardin(Q)
    return SDContext(
        P,
        Q,
        SweepStore(n_p, n_q),
        PrevMap(n_p, n_q),
        PriorityQueue{Int64,Float64}(),
        n_p,
        n_q,
    )
end


function get_event_points(c::SDContext{N,T,PT}, id::Int64) where {N,T,PT<:AbstractVector{T}}
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



# `prev_pts`: the points of the event _prev_id (get_event_points).
function SW_get_event_value(
    _id::Int64,
    prev_pts,
    prev_val::Float64,
    c::SDContext{N,T,PT},
) where {N,T,PT<:AbstractVector{T}}
    o_p, o_q = prev_pts
    n_p, n_q = get_event_points(c, _id)

    return prev_val + fp_match_price(o_p, n_p, o_q, n_q)
end


function is_sw_schedule_event(
    c::SDContext{N,T,PT},
    id::Int64,
)::Bool where {N,T,PT<:AbstractVector{T}}
    n_p, n_q = c.n_p, c.n_q
    if (EID_i(id) > n_p) || (EID_j(id) > n_q)
        return false
    end
    if (EID_i(id) >= n_p) && (! EID_i_is_vert(id))
        return false
    end
    if (EID_j(id) >= n_q) && (! EID_j_is_vert(id))
        return false
    end
    return ! sw_handled(c.store, id)
end


# Schedules event id, reached from id_prev; `prev_pts` = get_event_points(c, id_prev).
function SW_schedule_event!(
    c::SDContext{N,T,PT},
    id::Int64,
    (id_prev, prev_pts),
) where {N,T,PT<:AbstractVector{T}}
    if ! is_sw_schedule_event(c, id)
        return
    end
    st = c.store
    new_val = SW_get_event_value(id, prev_pts, sw_val(st, id_prev), c)
    if sw_val_or_inf(st, id) < new_val
        return  # A better event is already scheduled.
    end
    sw_set!(st, id, new_val, id_prev)

    if haskey(c.pq, id)
        delete!(c.pq, id)
    end
    push!(c.pq, id => new_val)
    return nothing
end



"""
    SweepDist_compute( P, Q )

A morphing between `P` and `Q` minimizing the approximate sweep price (the sum of
`fp_match_price` over matched edges), computed by a Dijkstra search over the
vertex-edge events. The minimum price is returned in `m.sol_value`; `m.leash` is
the longest leash of that morphing. The morphing need not be monotone.
"""
function SweepDist_compute(
    P::GPntSeq{N,T,PT},
    Q::GPntSeq{N,T,PT},
) where {N,T,PT<:AbstractVector{T}}
    check_curves(P, Q, "SweepDist_compute"; min_vertices = 2)
    c::SDContext{N,T,PT} = SDContext(P, Q)

    start_id = EID(1, true, 1, true)

    id_end = EID(c.n_p, true, c.n_q, true)
    st = c.store
    sw_set!(st, start_id, 0.0, start_id)
    push!(c.pq, start_id => 0.0)

    end_val = 0.0
    iters = 0
    pq = c.pq
    while ! isempty(pq)
        id = first(popfirst!(pq))

        if sw_handled(st, id)
            continue
        end

        sw_handle!(st, id)
        c.prev_map[id] = sw_prev(st, id)
        iters = iters + 1

        i = EID_i(id)
        j = EID_j(id)
        from = (id, get_event_points(c, id))

        if id == start_id
            SW_schedule_event!(c, EID(1, false, 1, true), from)
            SW_schedule_event!(c, EID(1, true, 1, false), from)
            continue
        end

        # Is it the *final* event?
        if (i == c.n_p) && (j == c.n_q)
            end_val = sw_val(st, id)
            break
        end

        # Is it on the boundary of the final cell?
        if is_final_cell(id, c.n_p, c.n_q)
            SW_schedule_event!(c, id_end, from)
            continue
        end
        SW_schedule_event!(c, EID(i+1, true, j, false), from)
        SW_schedule_event!(c, EID(i, false, j+1, true), from)
    end

    pes, qes = f_r_extract_solution_basic(P, Q, id_end, c.prev_map)

    morph::Morphing{N,T,PT} = Morphing_init(P, Q, pes, qes)
    morph.iters = iters
    morph.sol_value = end_val

    Morphing_verify_valid(morph)

    return morph
end


"""
    SweepDist_compute_refine_mono( P, Q )

A monotone morphing for the sweep distance. Runs `SweepDist_compute`; while the
result is not monotone, and monotonizing it would move the curves by 1% or more of
their total length, refines both curves along the morphing and runs it again.
Returns the monotone (or monotonized) morphing, whose curves can be refinements of
`P` and `Q`. Throws an error after `MAX_REFINEMENT_ROUNDS` rounds.
"""
function SweepDist_compute_refine_mono(
    poly_a::GPntSeq{N,T,PT},
    poly_b::GPntSeq{N,T,PT},
) where {N,T,PT<:AbstractVector{T}}
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

"""
    SweepDist_compute_split!( out, P, Q, iters, max_vertices = 1000000 )

Appends to `out` sweep-distance morphings at increasing resolution: first
`SweepDist_compute_refine_mono(P, Q)`, then `SweepDist_compute` on its curves with
the midpoint of every edge inserted (`split_edges`), repeated `iters` times or until
the two curves have more than `max_vertices` vertices together. Every morphing is
checked with `Morphing_verify_valid`. Returns `nothing`.
"""
function SweepDist_compute_split!(
    out::Vector{Morphing{N,T,PT}},
    poly_a::GPntSeq{N,T,PT},
    poly_b::GPntSeq{N,T,PT},
    iters::Int64,
    max_vertices::Int64 = 1000000,
) where {N,T,PT<:AbstractVector{T}}
    check_curves(poly_a, poly_b, "SweepDist_compute_split!"; min_vertices = 2)

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
and Q. The current implementation is naive, and very similar to the
discrete Frechet calculation.

"""
function SweepDist_lb_compute(
    P::GPntSeq{D,T,PT},
    Q::GPntSeq{D,T,PT},
) where {D,T,PT<:AbstractVector{T}}
    check_curves(P, Q, "SweepDist_lb_compute"; min_vertices = 2)

    pq = PriorityQueue{Tuple{Int64,Int64},Float64}()

    n_p = cardin(P)
    n_q = cardin(Q)

    # Inf marks a cell not reached yet (0 is a legitimate value).
    dp = fill(Inf, n_p, n_q)
    dp_dec_i = falses(n_p, n_q)
    handled = falses(n_p, n_q)

    # Every segment-segment distance is needed by up to two cells in each argument
    # order (the kernel is not symmetric bit for bit, so the orders are kept apart):
    # computed once, on first use. seg_pq[i, j]: P_i P_i+1 against Q_j Q_j+1; seg_qp[j, i]
    # the other way around. NaN: not computed yet.
    seg_pq = fill(NaN, n_p, n_q)
    seg_qp = fill(NaN, n_q, n_p)
    len_p = [fp_dist(P[i], P[i+1]) for i in 1:(n_p-1)]
    len_q = [fp_dist(Q[j], Q[j+1]) for j in 1:(n_q-1)]

    function seg_dist(memo, PA, i, PB, j)
        d = memo[i, j]
        isnan(d) || return d
        d = fp_dist_segment_segment(PA[i], PA[i+1], PB[j], PB[j+1])
        memo[i, j] = d
        return d
    end

    # edge_lb: for the grid edge (i, j) -> (i+1, j), the minimum of the distance
    # between the edge PA[i] PA[i+1] and the edges PB[j-1] PB[j] and PB[j] PB[j+1]
    # of the two adjacent cells: a lower bound on the leash along that edge.
    function edge_lb(memo, PA, i, PB, j)
        lb::Float64 = typemax(Float64)
        if i == cardin(PA)
            return 0.0
        end
        if j > 1
            lb = min(lb, seg_dist(memo, PA, i, PB, j - 1))
        end
        if j < cardin(PB)
            lb = min(lb, seg_dist(memo, PA, i, PB, j))
        end
        return lb
    end

    # The price of the grid edge (i, j) -> (i+1, j): the length of PA[i] PA[i+1] times
    # its edge_lb.
    function edge_v(memo, lens, (PA, i), (PB, j))::Float64
        (i == cardin(PA)) && return 0.0
        return lens[i] * edge_lb(memo, PA, i, PB, j)
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
    push!(pq, (1, 1) => dp[1, 1])

    iters::Int64 = 0
    while ! isempty(pq)
        iters = iters + 1

        ele = first(pq)  # the minimum (DataStructures 0.19)
        i::Int64, j::Int64 = ele[1]
        value = ele[2]

        popfirst!(pq)
        if handled[i, j]
            continue
        end

        handled[i, j] = true
        if ((i == n_p) && (j == n_q))
            break
        end


        # Now we need to schedule the next two adjacent entries..
        if (i < n_p)
            push_value(i, j, i + 1, j, value + edge_v(seg_pq, len_p, (P, i), (Q, j)))
        end
        if (j < n_q)
            push_value(i, j, i, j+1, value + edge_v(seg_qp, len_q, (Q, j), (P, i)))
        end
    end

    m = d_frechet_extract_solution(P, Q, dp_dec_i, n_p, n_q)
    m.iters = iters

    m.sol_value = dp[n_p, n_q]

    return m
end

# SweepDist.jl end of file
#############################################################################
