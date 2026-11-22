# Originally contributed by S. Har-Peled
# under MIT License

#---------------------------------------------------------
# Morphing
#
# A morphing is the matching between two polygons as computed by
# the Frechet distance.
#---------------------------------------------------------

include("DistFunc.jl");



"""
    FPointType

The kind of an `EventPoint`: `PT_VERTEX` (at a vertex of the curve) or `PT_ON_EDGE`
(inside an edge).
"""
@enum FPointType begin
    PT_VERTEX = 1
    PT_ON_EDGE = 2
end

"""
    EventPoint{N,T}

One end of a step of a `Morphing`: a location on a curve. `p` is the point; `i` is
the vertex (`type == PT_VERTEX`) or the first vertex of the edge (`type ==
PT_ON_EDGE`); `t` is the position along that edge, from 0 at vertex `i` to 1 at
vertex `i + 1`; `penalty` is the vertex offset charged to the event by computations
with offsets (`frechet_ve_r_compute_ext`), 0 otherwise. `EventPoint(p, i, type, t)`
sets the penalty to 0.
"""
struct EventPoint{N,T}
    p::Point{N,T}
    i::Int64
    type::FPointType
    t::Float64
    penalty::Float64
end

function EventPoint(_p::Point{N,T}, _i::Int64, _type::FPointType, _t::Float64) where {N,T}
    return EventPoint(_p, _i, _type, _t, 0.0)
end

function ev_update_p_t(ev::EventPoint{N,T}, p::Point{N,T}, t::Float64) where {N,T}
    return EventPoint(p, ev.i, ev.type, t, ev.penalty)
end
function ev_update_penalty(ev::EventPoint{N,T}, penalty) where {N,T}
    return EventPoint(ev.p, ev.i, ev.type, ev.t, penalty)
end
# leash = maximum length of an edge in the encoded matching.
#         (That is, the Frechet distance.)
"""
    Morphing

Encoding of a morphing (i.e., matching) between two polygonal curves.

`P` and `Q` are the curves the morphing was computed for (for `frechet_c_compute` and
`frechet_c_approx`, their refinements). They are shared with the inputs, not copied:
copy them before modifying them. `pes` and `qes` list the matched event points on `P`
and `Q`, and `leash` is the largest distance between matched points. Fields that an
algorithm does not compute stay at their defaults (`ratio`, `lower_bound`,
`monotone_err`, `sol_value`).

"""
@with_kw mutable struct Morphing{N,T}
    P::PntSeq{N,T}
    Q::PntSeq{N,T}
    pes::Vector{EventPoint{N,T}} # P event sequence
    qes::Vector{EventPoint{N,T}} # Q event sequence
    leash::Float64
    #    leash_offsets::Float64
    iters::Int64
    ratio::Float64
    monotone_err::Float64
    sol_value::Float64

    f_is_monotone_init::Bool
    f_is_monotone::Bool

    lower_bound::T
end

function Morphing_init(
    P::PntSeq{N,T},
    Q::PntSeq{N,T},
    pes::Vector{EventPoint{N,T}},
    qes::Vector{EventPoint{N,T}},
) where {N,T}

    @assert(length(pes) == length(qes))
    m = Morphing(
        P,
        Q,
        pes,
        qes,
        0.0,
        0,
        0.0,
        0.0,
        0.0,
        false,
        false,
        0.0, # lower_bound
    )
    Morphing_recompute_leash!(m)

    return m
end


"""
    Morphing_get_max_edges_err( m )

The pair of edges matched by `m` where the matching is loosest. For each step on
edges `P[i] P[i+1]` and `Q[j] Q[j+1]`, the error is the largest distance between
their endpoints minus the distance between the two segments. Returns `(i, j, err)`
for the largest error, or `(-1, -1, 0.0)` when no step has a positive one.
"""
function Morphing_get_max_edges_err(m::Morphing{D,T}) where {D,T}
    mx::Float64 = 0.0

    io = jo = -1

    P = m.P
    Q = m.Q

    l_p = cardin(P)
    l_q = cardin(Q)

    len = length(m.pes)
    for ind in 1:len
        i = m.pes[ind].i
        j = m.qes[ind].i

        if (i == l_p) || (j == l_q)
            continue
        end

        d = max(
            fp_dist(P[i], Q[j]),
            fp_dist(P[i], Q[j+1]),
            fp_dist(P[i+1], Q[j]),
            fp_dist(P[i+1], Q[j+1]),
        )
        ld = fp_dist_segment_segment(P[i], P[i+1], Q[j], Q[j+1])
        delta = d - ld
        if (delta > mx)
            mx = delta
            io, jo = i, j
        end
    end

    return io, jo, mx
end


####################################################
# Compute for every vertex the maximum leash length used with this
# vertex.
####################################################

function extract_vertex_leash(
    P::PntSeq{N,T},
    pes::Vector{EventPoint{N,T}},
    qes::Vector{EventPoint{N,T}},
) where {N,T}
    P_len = cardin(P)
    len = length(pes)
    vl = zeros(Float64, P_len)
    for i in 1:len
        ev = pes[i]
        evx = qes[i]
        # BUG: The if should not have been there... We are being very
        # conservative about vertices that should not be simplified.

        # if ( ev.type == PT_VERTEX )
        vl[ev.i] = max(vl[ev.i], fp_dist(ev.p, evx.p))
        #end
    end

    return vl
end


"""
    Morphing_extract_vertex_radii

    Computes for each polygon vertex, the length of the longest edge
    in the matching attached to it. It is a cheap upper bound on the
    local Frechet distance for each vertex (and implicitly the
    attached edge).
"""
function Morphing_extract_vertex_radii(m::Morphing{N,T}) where {N,T}
    pl = extract_vertex_leash(m.P, m.pes, m.qes)
    ql = extract_vertex_leash(m.Q, m.qes, m.pes)
    return pl, ql
end


function Morphing_swap_sides!(m::Morphing{N,T}) where {N,T}
    m.P, m.Q = m.Q, m.P
    #    swap!( m.P, m.Q );
    #swap!( m.pes, m.qes );
    m.pes, m.qes = m.qes, m.pes
end

"""
    Morphing_empty

Returns an empty morphing (i.e., a constructor).
"""
function Morphing_empty(P::PntSeq{N,T}, Q::PntSeq{N,T}) where {N,T}
    check_coordinate_type(P, "Morphing_empty", "P")
    check_coordinate_type(Q, "Morphing_empty", "Q")
    pes = Vector{EventPoint{N,T}}()
    qes = Vector{EventPoint{N,T}}()
    r::Float64 = -1
    return Morphing(P, Q, pes, qes, r, 0, 0.0, 0.0, 0.0, false, false, 0.0)
end


###################################################################
# Check if the matching is monotone. That is, the endpoints of the
# matching on each curve are sorted, in the same order as the
# matching.
##################################################################
function events_seq_is_monotone(s::Vector{EventPoint{N,T}}) where {N,T}

    len = length(s)
    i::Int64 = 1
    while (i <= len)
        ep = s[i]
        if ep.type == PT_VERTEX
            i = i + 1
            continue
        end

        j = i
        loc = s[i].i
        while ((j < len) && (s[j+1].type == PT_ON_EDGE) && (s[j+1].i == loc))
            j = j + 1
        end

        # i:j is the sequence of edge-vertex events
        for k in i:(j-1)
            if s[k+1].t < s[k].t
                return false
            end
        end
        i = j + 1
    end

    return true
end


##############################################################################
# Make the matching (which might be potentially not monotone) into a monotone
# one.
##############################################################################
function events_seq_make_monotone(P::PntSeq{N,T}, s::Vector{EventPoint{N,T}}) where {N,T}
    ns = Vector{EventPoint{N,T}}()

    len = length(s)
    i::Int64 = 1
    delta = 0.0
    while (i <= len)
        ep = s[i]
        if ep.type == PT_VERTEX
            up = deepcopy(ep)
            push!(ns, up)
            i = i + 1
            continue
        end

        j = i
        loc = s[i].i
        t = s[i].t
        while ((j < len) && (s[j+1].type == PT_ON_EDGE) && (s[j+1].i == loc))
            nep = deepcopy(s[j])
            t = max(t, nep.t)
            if (t > nep.t)
                ell = fp_dist(P[nep.i], P[nep.i+1])
                delta = max(delta, (t - nep.t) * ell)
                p = fp_convex_comb(P[nep.i], P[nep.i+1], t)
                nep = ev_update_p_t(nep, p, t)
            end
            push!(ns, nep)
            j = j + 1
        end

        nep = deepcopy(s[j])
        if (t > nep.t)
            ell = fp_dist(P[nep.i], P[nep.i+1])
            delta = max(delta, (t - nep.t) * ell)
            p = fp_convex_comb(P[nep.i], P[nep.i+1], t)
            nep = ev_update_p_t(nep, p, t)
        end
        push!(ns, nep)

        i = j + 1
    end

    @assert(length(ns) == length(s))
    return ns, delta
end


##############################################################################
# The leash of the matching after it is made monotone (as by
# events_seq_make_monotone), without building the monotone matching.
##############################################################################
function events_seq_get_monotone_leash(
    P::PntSeq{N,T},
    s::Vector{EventPoint{N,T}},
    s_alt::Vector{EventPoint{N,T}},
    leash::T,
) where {N,T}

    len = length(s)
    i::Int64 = 1
    while (i <= len)
        ep = s[i]
        if ep.type == PT_VERTEX
            i = i + 1
            continue
        end

        j = i
        loc = s[i].i
        t = s[i].t
        while ((j < len) && (s[j+1].type == PT_ON_EDGE) && (s[j+1].i == loc))
            nep = s[j]
            t = max(t, nep.t)
            if (t > nep.t)
                ell = fp_dist(P[nep.i], P[nep.i+1])
                new_p = fp_convex_comb(P[nep.i], P[nep.i+1], t)
                leash = max(leash, fp_dist(new_p, s_alt[j].p))
            end
            j = j + 1
        end

        nep = s[j]
        if (t > nep.t)
            ell = fp_dist(P[nep.i], P[nep.i+1])
            new_p = fp_convex_comb(P[nep.i], P[nep.i+1], t)
            leash = max(leash, fp_dist(new_p, s_alt[j].p))
        end
        i = j + 1
    end

    return leash
end


function Morphing_recompute_leash!(m::Morphing{N,T}) where {N,T}
    @assert(length(m.pes) == length(m.qes))
    r = 0.0
    for i in eachindex(m.pes)
        ell = (fp_dist(m.pes[i].p, m.qes[i].p) - m.pes[i].penalty)
        if ell > r
            r = ell
        end
    end
    # exact: the same operations recompute the same leash; 0 means not computed yet.
    if (m.leash != r) && (m.leash != 0)
        error(
            "Morphing_recompute_leash: stored leash $(m.leash) differs from the recomputed leash $r",
        )
    end
    m.leash = r
end

"""
    Morphing_as_polygons

Turns the morphing matching into a "true" matching, by creating two
polygons whose edges are directly matched.  The output polygons P
and Q will definitely have repeated points.

"""
function Morphing_as_polygons(m::Morphing{N,T}) where {N,T}
    P = PntSeq{N,T}()
    Q = PntSeq{N,T}()

    Morphing_verify_valid(m)

    #leash::Float64 = 0.0;
    for i in eachindex(m.pes)
        push!(P, m.pes[i].p)
        push!(Q, m.qes[i].p)
    end

    return P, Q
end


"""
    Morphing_as_polygons_w_times( m )

The morphing as two polygons with one vertex per step (`Morphing_as_polygons`),
and a time in [0, 1] for each step: the distance the two points have travelled
together so far, divided by the total. Returns `(P, Q, t)`. A step where neither
point moves takes no time, so times can repeat. If nothing moves at all, the times
are spread evenly.
"""
function Morphing_as_polygons_w_times(m::Morphing{N,T}) where {N,T}
    P, Q = Morphing_as_polygons(m)
    travelled = zeros(length(P))
    for i in 1:(length(P)-1)
        travelled[i+1] = travelled[i] + fp_dist(P[i], P[i+1]) + fp_dist(Q[i], Q[i+1])
    end
    total = last(travelled)
    if total == 0.0  # exact: nothing moves
        # Every time in [0, 1] must fall between two events.
        return P, Q, collect(range(0.0, 1.0, length = length(P)))
    end
    # Nondecreasing, and the last time is total / total == 1.0 exactly.
    return P, Q, travelled ./ total
end


"""
    Morphing_as_function_w_times( m )

Like `Morphing_as_polygons_w_times`, but the times follow the arc length of `P`
alone. Returns `(P, Q, t)`: the morphing as two polygons with one vertex per step,
and a time in [0, 1] for each vertex. Steps where `P` does not move get a small
positive duration, so `t` increases strictly; a constant `P` gets all times 0.
"""
function Morphing_as_function_w_times(m::Morphing{N,T}) where {N,T}
    P, Q = Morphing_as_polygons(m)
    t = Vector{Float64}()

    len = geom_length(P)
    if len == 0.0  # exact: a constant curve; times are relative to its length
        return P, Q, zeros(length(P))
    end
    push!(t, 0.0)
    c_z = 0
    fake_delta = len / (4 + length(P))
    for i in 1:(length(P)-1)
        delta = fp_dist(P[i], P[i+1]) / len
        if (delta == 0.0)  # exact: a step where P does not move
            c_z = c_z + 1
            len = len + fake_delta
        end
    end
    for i in 1:(length(P)-1)
        if (fp_dist(P[i], P[i+1]) == 0.0)  # exact: a step that does not move
            delta = fake_delta / len
        else
            delta = fp_dist(P[i], P[i+1]) / len
        end

        push!(t, min(last(t) + delta, 1.0))
    end
    if last(t) != 1.0  # exact: rounding can leave the last time just below 1
        pop!(t)
        push!(t, 1.0)
    end
    return P, Q, t
end


"""
    polygons_get_loc_at_time( P, Q, times, t )

The two locations at time `t` of a morphing given as polygons with one vertex per
step and their `times` (as from `Morphing_as_polygons_w_times`): interpolates
between the two vertices whose times bracket `t`. Returns `(p, q)`. Times outside
[0, 1] give the first or the last pair.
"""
function polygons_get_loc_at_time(
    P::PntSeq{D,T},
    Q::PntSeq{D,T},
    times::Vector{T},
    t::Real,
) where {D,T}
    (t <= 0.0) && return first(P), first(Q)
    (t >= 1.0) && return last(P), last(Q)

    pos = searchsortedfirst(times, t)
    (1 < pos <= length(P)) ||
        error("polygons_get_loc_at_time: time $t is outside the event times")

    # searchsortedfirst gives times[pos-1] < t <= times[pos], so the denominator is
    # positive. Repeated times (zero-duration steps of the morphing) need no special
    # case: the bracket never straddles them.
    prev = pos - 1
    delta = (t - times[prev]) / (times[pos] - times[prev])
    p = fp_convex_comb(P[prev], P[pos], delta)
    q = fp_convex_comb(Q[prev], Q[pos], delta)
    return p, q
end

"""
    Morphing_sample_uniformly( m, n )

Samples the morphing at `n >= 2` evenly spaced times in [0, 1] (the times of
`Morphing_as_polygons_w_times`). Returns an `n × 2N` matrix: row `k` holds the point
on `P` (columns `1:N`) and the matching point on `Q` (columns `N+1:2N`).
"""
function Morphing_sample_uniformly(m::Morphing{N,T}, n::Int64) where {N,T}
    n >= 2 ||
        throw(ArgumentError("Morphing_sample_uniformly: need at least 2 samples (got $n)"))
    P, Q, times = Morphing_as_polygons_w_times(m)

    out = zeros(n, 2*N)
    for i in 0:(n-1)
        t = i/(n-1)
        p, q = polygons_get_loc_at_time(P, Q, times, t)
        out[i+1, 1:N] .= p
        out[i+1, (N+1):(2N)] .= q
    end

    return out
end

function Morphing_is_monotone(m::Morphing{N,T}) where {N,T}
    if m.f_is_monotone_init
        return m.f_is_monotone
    end
    m.f_is_monotone = (events_seq_is_monotone(m.pes) && events_seq_is_monotone(m.qes))
    m.f_is_monotone_init = true

    return m.f_is_monotone
end

##########################################################3
# Turn the morphing into a monotone morphing...
##########################################################3
"""
    Morphing_monotonize(m)

Turns a morphing into a monotone morphing, by simply not going back,
staying in place if necessary. `monotone_err` of the result is how far a point
had to be held back.

The result may share memory with `m`: if `m` is already monotone, it is returned
itself, not a copy; otherwise the result refers to the same curves `m.P` and
`m.Q`. Use `deepcopy` before changing either.
"""
function Morphing_monotonize(m::Morphing{N,T}) where {N,T}
    if Morphing_is_monotone(m)
        return m  # aliased, by design: see the docstring
    end
    P = m.P
    Q = m.Q

    pes_new, da = events_seq_make_monotone(P, m.pes)
    qes_new, db = events_seq_make_monotone(Q, m.qes)


    m_out = Morphing_init(P, Q, pes_new, qes_new)

    m_out.monotone_err = max(da, db)

    return m_out
end


"""
    Morphing_monotone_leash( m )

The leash of `Morphing_monotonize(m)`, computed without building the monotone
morphing; `m.leash` when `m` is already monotone.
"""
function Morphing_monotone_leash(m::Morphing{N,T}) where {N,T}
    if Morphing_is_monotone(m)
        return m.leash # deepcopy( m );
    end

    leash = m.leash
    leash = events_seq_get_monotone_leash(m.P, m.pes, m.qes, leash)
    leash = events_seq_get_monotone_leash(m.Q, m.qes, m.pes, leash)

    return leash
end



##########################################################################


function check_times(V::Vector{EventPoint{N,T}}) where {N,T}
    for ev in V
        if (ev.t < 0.0) || (ev.t > 1.0)
            throw(ArgumentError("morphing event time $(ev.t) is outside [0, 1]"))
        end
    end
end

function check_no_nan(V::Vector{EventPoint{N,T}}) where {N,T}
    for ev in V
        any(isnan, ev.p) &&
            throw(ArgumentError("morphing event has a NaN coordinate: $(ev.p)"))
    end
end


"Whether `p` and `q` agree up to a rounding tolerance relative to the coordinates."
function event_close(p::Point{N,T}, q::Point{N,T}) where {N,T}
    scale = max(maximum(abs, p), maximum(abs, q))
    return fp_dist(p, q) <= 1.0e-9 * (1.0 + scale)
end

"""
    check_events( C, V, name )

Throws an `ArgumentError` unless the event sequence `V` runs from the first to the
last vertex of the curve `C`, every index is a vertex (or, for an edge event, the
start of an edge) of `C`, and every location agrees with its index and coefficient.
"""
function check_events(C::PntSeq{N,T}, V::Vector{EventPoint{N,T}}, name) where {N,T}
    n = cardin(C)
    isempty(V) && throw(ArgumentError("morphing: the $name event sequence is empty"))
    (event_close(first(V).p, C[1]) && event_close(last(V).p, C[n])) || throw(
        ArgumentError(
            "morphing: the $name events do not run from the first to the last vertex",
        ),
    )
    for ev in V
        last_ok = (ev.type == PT_VERTEX) ? n : n - 1
        (1 <= ev.i <= last_ok) || throw(
            ArgumentError("morphing: $name event index $(ev.i) is outside 1:$last_ok"),
        )
        loc = (ev.type == PT_VERTEX) ? C[ev.i] : fp_convex_comb(C[ev.i], C[ev.i+1], ev.t)
        event_close(ev.p, loc) || throw(
            ArgumentError(
                "morphing: $name event at $(ev.p) does not lie at index $(ev.i), t = $(ev.t)",
            ),
        )
    end
    return nothing
end

"""
    Morphing_verify_valid( m )

Throws an `ArgumentError` unless `m` is structurally valid: `pes` and `qes` have the
same length; each runs from the first to the last vertex of its curve; every event has
a coefficient in [0, 1], no NaN, an index within its curve, and a location that agrees
with its index and coefficient. Monotonicity is not required (see
`Morphing_is_monotone`).
"""
function Morphing_verify_valid(m::Morphing{N,T}) where {N,T}
    length(m.pes) == length(m.qes) || throw(
        ArgumentError("morphing: $(length(m.pes)) events on P but $(length(m.qes)) on Q"),
    )
    check_times(m.pes)
    check_times(m.qes)

    check_no_nan(m.pes)
    check_no_nan(m.qes)

    check_events(m.P, m.pes, "P")
    check_events(m.Q, m.qes, "Q")
    return nothing
end


###########################################################################
# Compute the morphing src(u) -> trg(v): that is, v(u(t))
###########################################################################

"""
    Morphing_combine( u, v )

Takes two morphings u, v (i.e., two parameterizations) and combines them
into a single morphing v(u(⋅)): first u, then v.

For example, if u: γ → δ  and  v: δ → ξ, then the returned morphing is
v(u(⋅)): γ → ξ, from `u.P` to `v.Q`. The middle curves `u.Q` and `v.P` must
have the same length (they are normally the same curve).

"""
function Morphing_combine(u::Morphing{N,T}, v::Morphing{N,T}) where {N,T}
    u_prm = Morphing_extract_prm(u)
    v_prm = Morphing_extract_prm(v)

    (cardin(u_prm) > 1 && cardin(v_prm) > 1) ||
        throw(ArgumentError("Morphing_combine: both morphings need at least two events"))

    prm = parameterization_combine(v_prm, u_prm)

    pes, qes = event_sequences_extract(prm, u.P, v.Q)

    m = Morphing_init(u.P, v.Q, pes, qes)

    return m
end

"""
    Morphing_SweepDist_approx_price( m )

A fast estimate of `Morphing_SweepDist_price(m)`: for every step of the morphing,
the two segment lengths times the average of the leash lengths at the step's ends
(`fp_match_price`), summed.
"""
function Morphing_SweepDist_approx_price(m::Morphing{N,T}) where {N,T}
    P, Q = Morphing_as_polygons(m)

    @assert(length(P) == length(Q))
    len = cardin(P)
    price::Float64 = 0
    for i in 1:(len-1)
        price = price + fp_match_price(P[i], P[i+1], Q[i], Q[i+1])
    end

    return price
end


"""
    Morphing_SweepDist_price( m )

The sweep price of the morphing `m`: the sum, over its steps, of the sweep distance
between the two matched segments (`SweepDist_segs`). Throws an `ArgumentError` if
`m` is not valid (`Morphing_verify_valid`).
"""
function Morphing_SweepDist_price(m::Morphing{N,T}) where {N,T}
    Morphing_verify_valid(m)

    P, Q = Morphing_as_polygons(m)

    @assert(length(P) == length(Q))
    len = cardin(P)
    price::Float64 = 0.0
    for i in 1:(len-1)
        delta::Float64 = SweepDist_segs(P[i], P[i+1], Q[i], Q[i+1])
        price = price + delta
    end

    return price
end

"""
    morphing_profile( m, n )

The leash length along the morphing, at `n` samples. Returns `(X, Y)`: `X` runs
evenly from 0 to the total distance `L` that the two points travel (the lengths of
the two polygons of `Morphing_as_polygons`), and `Y[k]` is the leash length once the
points have travelled `X[k]` together: at the time `(k-1)/(n-1)` (`= X[k] / L`) of
`Morphing_as_polygons_w_times`.
"""
function morphing_profile(m::Morphing{D,T}, n::Int64) where {D,T}

    PB, QB, times = Morphing_as_polygons_w_times(m)
    # The times of Morphing_as_polygons_w_times measure the distance both points travel.
    len = geom_length(PB) + geom_length(QB)
    X = range(0, len, length = n)

    Y = Vector{Float64}()
    for k in 1:n
        # The time of X[k], computed directly: X[k] / len rounds, and is NaN for len = 0.
        pos = (n == 1) ? 0.0 : (k - 1) / (n - 1)
        p, q = polygons_get_loc_at_time(PB, QB, times, pos)
        push!(Y, fp_dist(p, q))
    end
    return X, Y
end


const Morphing2F = Morphing{2,Float64};

#
# End of file
##########################################################################
