###########################################################################
# Compute the morphing src(v) -> trg(u): That is u(v(t))
###########################################################################

"""
    frechet_dist_upper_bound

Returns a rough upper bound on the Frechet distance between the two
curves. This upper bound is on the continuous distance. There is no guarantee
on how bad the approximation is. This is used as a starting point for
real approximation of the Frechet distance, and should not be used
otherwise.

"""
function frechet_dist_upper_bound(poly_a::PntSeq{N,T}, poly_b::PntSeq{N,T}) where {N,T}
    check_curves(poly_a, poly_b, "frechet_dist_upper_bound"; min_vertices = 2)
    w_a = frechet_width_approx(poly_a)
    w_b = frechet_width_approx(poly_b)
    w = max(fp_dist(poly_a[1], poly_b[1]), fp_dist(last(poly_a), last(poly_b)))
    return w_a + w_b + w
end


"""
    frechet_ve_r_mono_compute

Computes the VE-retractable Frechet distance, and then monotonizes it. No
guarantee on the quality of the distance this realizes.

"""
function frechet_ve_r_mono_compute(poly_a::PntSeq{D,T}, poly_b::PntSeq{D,T}) where {D,T}
    m = frechet_ve_r_compute(poly_a, poly_b)
    mm = Morphing_monotonize(m)

    return mm
end


# Lifts the morphing m between the simplified curves back to a monotone morphing
# between the originals, given the simplification morphisms m_p (poly_a to P) and
# m_q (poly_b to Q). Also returns the total simplification error.
function combine_via_simplified(
    m_p::Morphing{N,T},
    m::Morphing{N,T},
    m_q::Morphing{N,T},
) where {N,T}
    Morphing_swap_sides!(m_q)

    # m = m_p * mm * m_q
    mmu = Morphing_monotonize(Morphing_combine(m_p, m))
    mm_out = Morphing_combine(mmu, m_q)

    return mm_out, m_p.leash + m_q.leash
end


# Records a certified lower bound `lb` on the Frechet distance in `m`, and the
# ratio between the leash (an upper bound) and it. `lb` is clamped to the leash, which
# it can exceed only by rounding.
function set_bounds!(m::Morphing, lb::Real)
    m.lower_bound = clamp(lb, 0.0, m.leash)
    m.ratio = (m.lower_bound > 0.0) ? m.leash / m.lower_bound : (m.leash == 0.0 ? 1.0 : Inf)  # exact: a zero leash is optimal
    return m
end


# Upper/lower ratio certified for a morphing built from simplifications with
# total error total_err.
function approx_ratio_bound(m_leash, t_approx, total_err, out_leash, approx)
    lbx = (m_leash / t_approx) - total_err
    if (total_err == 0.0) && (m_leash == 0.0)  # exact: identical curves, ratio 1
        return 1.0
    end
    if (lbx <= 0.0)
        return 2.0 + approx # Very bad.
    end
    return out_leash / lbx
end


"""
    frechet_c_approx

Approximates the continuous Frechet distance between the two input
curves. Returns a monotone morphing realizing it.

# Arguments

- `approx` : The output morphing has Frechet distance <= approx*optimal.

Importantly, approx can be larger than 2, if you want a rough
approximation.

The quality of approximation is available at ret.ratio: ret.lower_bound
(= ret.leash/ret.ratio) is a certified lower bound on the Frechet distance,
while ret.leash is an upper bound.

"""
function frechet_c_approx(
    poly_a::PntSeq{N,T},
    poly_b::PntSeq{N,T},
    approx::Real,
)::Morphing{N,T} where {N,T}
    check_curves(poly_a, poly_b, "frechet_c_approx"; min_vertices = 2)
    approx = check_approx(approx, "frechet_c_approx")


    t_approx::Float64 = approx

    f_do_one_round::Bool = true
    d = frechet_dist_upper_bound(poly_a, poly_b)


    # r: Radius of simplification allowed
    r::Float64 = d/(approx + 4.0) #max( d/4, d * min( approx, 1.0 ) );

    P::PntSeq{N,T} = poly_a
    Q::PntSeq{N,T} = poly_b

    p_indices = Vector{Int64}()
    q_indices = Vector{Int64}()

    m = Morphing_empty(poly_a, poly_b)
    rounds::Int = 0
    while (true)
        rounds += 1
        check_rounds(rounds, "frechet_c_approx")
        while ((r >= (d / (approx + 4.0))) || f_do_one_round)
            r = r / 2.0
            P, p_indices = Polygon_simplify_ext(poly_a, r)
            Q, q_indices = Polygon_simplify_ext(poly_b, r)

            if ((2 * cardin(P) > cardin(poly_a)) && (2 * cardin(Q) > cardin(poly_b)))
                # Its lower bound is the VE-retractable distance, not the leash.
                m = frechet_mono_via_refinement(poly_a, poly_b, min(approx, 1.00001))[1]
                return set_bounds!(m, m.lower_bound)
            end

            @assert((cardin(P) > 1) && (cardin(Q) > 1))

            t_approx = 3.0/4.0 + t_approx / 4.0

            m = frechet_mono_via_refinement(P, Q, t_approx)[1]
            d = m.leash  #        frechet_ve_r_compute_dist( P, Q )
            if (fp_equal(d, 0.0))
                break
            end

            # d - 2r: Lower bound on the real Frechet distance.
            f_do_one_round = false
        end

        mm_out, total_err = combine_via_simplified(
            frechet_c_mono_approx_subcurve(poly_a, P, p_indices)[1],
            m,
            frechet_c_mono_approx_subcurve(poly_b, Q, q_indices)[1],
        )
        ratio = approx_ratio_bound(m.leash, t_approx, total_err, mm_out.leash, approx)
        if ratio <= approx
            # m.leash <= t_approx d(P, Q) and d(P, Q) <= d(poly_a, poly_b) + total_err.
            return set_bounds!(mm_out, m.leash / t_approx - total_err)
        end
        f_do_one_round = true
    end
end




function Polygon_simplify_radii_ext(P::PntSeq{N,T}, r::Vector{T}) where {N,T}

    PS, p_indices = Polygon_simplify_radii(P, r)
    @assert(length(p_indices) == cardin(PS))
    if (cardin(PS) <= (cardin(P) / 2))
        return PS, p_indices, false
    end

    p_indices = Vector{Int64}()
    append!(p_indices, 1:cardin(P))
    return P, p_indices, true
end


function propagate_mins!(qz, rounds::Int64 = 1)
    a = deepcopy(qz)
    for r in 1:rounds
        for i in eachindex(a)
            if (i == 1) || (i == length(qz))
                continue
            end

            a[i] = min(qz[i-1], qz[i], qz[i+1])
        end

        for i in eachindex(a)
            qz[i] = a[i]
        end
    end
end

function eq(a::Float64, b::Float64, tolerance::Float64)::Bool
    if a == b
        return true
    end
    return abs(a - b) <= (tolerance * (abs(a) + abs(b)))
end


function get_refinement_map(PSR::PntSeq{N,T}, PS::PntSeq{N,T}) where {N,T}
    @assert(cardin(PSR) >= cardin(PS))

    len = cardin(PSR)
    lenPS = cardin(PS)
    mp = zeros(Int64, len)

    curr = 1
    i = 1
    while (curr <= len)
        if (i < lenPS) && (PSR[curr] == PS[i+1])
            i = i + 1
        end
        mp[curr] = i
        curr = curr + 1
    end

    return mp
end


function offsets_copy_map(map_PS, offsets)
    offs = zeros(Float64, length(map_PS))
    for i in eachindex(map_PS)
        offs[i] = offsets[map_PS[i]]
    end
    return offs
end


"""
    simplify_morphing_sensitive( m, factor )

Simplifies both curves of `m` where the morphing leaves slack: a vertex may move by
up to `(m.lower_bound - r) / factor`, where `r` is the leash length at the vertex;
vertices without slack, and their neighbours, stay. Returns the simplified curves
`(PS, QS)`. `factor` must be positive.
"""
function simplify_morphing_sensitive(m::Morphing{D,T}, factor::Float64) where {D,T}
    (isfinite(factor) && factor > 0.0) || throw(
        ArgumentError(
            "simplify_morphing_sensitive: factor must be finite and positive (got $factor)",
        ),
    )
    pl, ql = Morphing_extract_vertex_radii(m)

    lower_bound = m.lower_bound

    pz = ((lower_bound * ones(length(pl))) - pl) / factor
    qz = ((lower_bound * ones(length(ql))) - ql) / factor

    propagate_mins!(pz, 1)
    propagate_mins!(qz, 1)


    PS, p_indices, f_PS_exact = Polygon_simplify_radii_ext(m.P, pz)
    QS, q_indices, f_QS_exact = Polygon_simplify_radii_ext(m.Q, qz)

    return PS, QS
end


# Tuning parameters of frechet_c_compute.
const FC_ROUGH_APPROX = 2.0           # factor of the first, rough approximation
const FC_MIN_APPROX = 1.00000000001   # smallest approximation factor ever requested
const FC_SLACK_DIVISOR = 4.0          # vertex slack = (bound - vertex radius) / divisor
const FC_SLACK_GROWTH = 1.3           # after a failed round the divisor grows by this
const FC_APPROX_SHRINK = 10.0         # ... and (approx - 1) shrinks by this
const FC_SLACK_ROUNDS = 8             # smoothing rounds of propagate_mins!

# How far each vertex may move in the simplification: (bound - radius) / divisor,
# where `radii` are the leash lengths at the vertices. Smoothed so that the
# neighbours of a vertex that must not move (negative slack) do not move either.
function vertex_slack(radii::Vector{Float64}, bound::Float64, divisor::Float64)
    slack = ((bound * ones(length(radii))) - radii) / divisor
    propagate_mins!(slack, FC_SLACK_ROUNDS)
    return slack
end

# Turns a morphing `m_mid` between the simplifications PS and QS into a morphing
# between P and Q, through P -> PS -> QS -> Q. Each side is (curve, simplification,
# indices of the simplification's vertices in the curve). Also returns the two
# simplification morphings and their offsets.
function stitch_simplified(
    m_mid::Morphing{N,T},
    p_side::Tuple{PntSeq{N,T},PntSeq{N,T},Vector{Int64}},
    q_side::Tuple{PntSeq{N,T},PntSeq{N,T},Vector{Int64}},
) where {N,T}
    P, PS, p_indices = p_side
    Q, QS, q_indices = q_side
    m_a, offsets_a = frechet_c_mono_approx_subcurve(P, PS, p_indices)
    m_pq = Morphing_combine(m_a, m_mid)
    m_b, offsets_b = frechet_c_mono_approx_subcurve(Q, QS, q_indices)
    Morphing_swap_sides!(m_b)
    return Morphing_combine(m_pq, m_b), m_a, m_b, offsets_a, offsets_b
end

# The VE-retractable distance with offsets between the refined simplifications PSR
# and QSR, a lower bound on d(P, Q). Each side is (refined simplification,
# simplification, offsets at the simplification's vertices); the offsets are carried
# over to the refinement. The refined curves must be used here, so that the distance
# is a monotone one.
function offsets_lower_bound(
    p_side::Tuple{PntSeq{N,T},PntSeq{N,T},Vector{Float64}},
    q_side::Tuple{PntSeq{N,T},PntSeq{N,T},Vector{Float64}},
) where {N,T}
    PSR, PS, offsets_a = p_side
    QSR, QS, offsets_b = q_side
    PS_offs = offsets_copy_map(get_refinement_map(PSR, PS), offsets_a)
    QS_offs = offsets_copy_map(get_refinement_map(QSR, QS), offsets_b)
    @assert(length(PS_offs) == cardin(PSR))
    @assert(length(QS_offs) == cardin(QSR))
    return frechet_ve_r_compute_ext(PSR, QSR, PS_offs, QS_offs, true)
end

###########################################################################
###########################################################################
"""
    frechet_c_compute

Compute the exact continuous (monotone) Frechet distance between the
two polygons. It should be reasonably fast.

This function is somewhat slower than the approximate versions. Use it
only if you really want the exact answer. Consider using
`frechet_c_approx` instead.

The result also carries a certified lower bound, `m.lower_bound`, and
`m.ratio = m.leash / m.lower_bound` (1 when the distance is exact).

# Details

This works by first computing a very rough approximation, followed by
distance-sensitive simplification of the curves. It then computes the
monotone fr_ve_r distance between the simplified curves, and
combines it to get a distance between the two original curves. It makes
sure the answers are the same, otherwise, it repeats with a finer
simplification/approximation until they are equal.

Finally, the algorithm uses the fr_ve_r_with_offsets distance between
the two simplified curves to compute a lower bound, and makes sure this
is equal to the Frechet distance computed. If they are equal, then the
upper/lower bounds on the Frechet distance of the two curves are the
same, which implies that the computed distance is indeed the desired
Frechet distance.

# More details

To ensure convergence, the monotone distance between the simplified curves is
computed by refinement (`frechet_mono_via_refinement`), so that it approaches their
VE-retractable distance.

# Arguments
- f_accept_approx: return once the leash is within a factor `1 + tolerance` of
  the certified lower bound. The test for agreement up to `tolerance` runs first
  and accepts slightly more, so this flag matters only when rounding puts the
  lower bound above the leash. It is kept for compatibility.

- Tolerance: If |a-b| <= tolerance*( |a| + |b| ) then we consider
  numbers to be equal. In a perfect world tolerance would be zero.
  But floating point issues require us to pick some value.

"""
function frechet_c_compute(
    P::PntSeq{N,T},
    Q::PntSeq{N,T},
    f_accept_approx::Bool = true,
    tolerance::Real = 0.00001,
) where {N,T}
    check_curves(P, Q, "frechet_c_compute"; min_vertices = 2)
    tolerance = check_tolerance(tolerance, "frechet_c_compute")

    # Repeated consecutive vertices describe the same curve. The refinement below
    # drops them, which breaks the vertex correspondence it relies on
    # (get_refinement_map), so they are removed up front.
    P = remove_consecutive_duplicates(P)
    Q = remove_consecutive_duplicates(Q)

    aprx_refinement::Float64 = FC_ROUGH_APPROX
    divisor::Float64 = FC_SLACK_DIVISOR

    m_rough = frechet_c_approx(P, Q, aprx_refinement)
    if (m_rough.leash == 0.0)  # exact: the curves coincide; nothing to refine
        return m_rough
    end

    # A good approximation, at most 1.1, finer when the curves are long compared to
    # their distance. Curves far apart relative to their length round it to 1.
    approx = min(1.0 + (geom_length(P) + geom_length(Q)) / (100.0 * m_rough.leash), 1.1)
    if approx <= 1.0
        approx = FC_MIN_APPROX
    end
    m = (m_rough.ratio > approx) ? frechet_c_approx(P, Q, approx) : m_rough

    # Leash lengths at the vertices, which decide how far each may move.
    pl, ql = Morphing_extract_vertex_radii(m)
    # Only steers the simplification; certified_lb is the bound that is reported.
    slack_bound = m.leash / m.ratio
    certified_lb = m.lower_bound

    round::Int64 = 0
    while true
        round = round + 1
        check_rounds(round, "frechet_c_compute")
        aprx_refinement = max(aprx_refinement, FC_MIN_APPROX)

        # Simplify both curves, moving each vertex at most by its slack.
        PS, p_indices, _ =
            Polygon_simplify_radii_ext(P, vertex_slack(pl, slack_bound, divisor))
        QS, q_indices, _ =
            Polygon_simplify_radii_ext(Q, vertex_slack(ql, slack_bound, divisor))

        # Simplifications this large save nothing: compute the distance directly.
        if ((2 * cardin(PS) > cardin(P)) && (2 * cardin(QS) > cardin(Q)))
            m = frechet_mono_via_refinement(P, Q, 1.0 + tolerance)[1]
            return set_bounds!(m, max(certified_lb, m.lower_bound))
        end

        t_aprx = max(tolerance + 1.0, aprx_refinement)
        m_mid, _, PSR, QSR = frechet_mono_via_refinement(PS, QS, t_aprx)
        mw, m_a, m_b, offsets_a, offsets_b =
            stitch_simplified(m_mid, (P, PS, p_indices), (Q, QS, q_indices))
        # m_mid.lower_bound <= d(PS, QS) <= d(P, Q) + m_a.leash + m_b.leash.
        certified_lb = max(certified_lb, m_mid.lower_bound - m_a.leash - m_b.leash)

        # The simplification increased the leash: simplify less and retry.
        if (!eq(m_mid.leash, mw.leash, tolerance))
            slack_bound = max(slack_bound, mw.leash - 2.0*m_a.leash - 2.0*m_b.leash)
            divisor = divisor * FC_SLACK_GROWTH
            aprx_refinement = 1.0 + (aprx_refinement - 1.0) / FC_APPROX_SHRINK
            pl, ql = Morphing_extract_vertex_radii(mw)
            continue
        end

        m_final = offsets_lower_bound((PSR, PS, offsets_a), (QSR, QS, offsets_b))
        certified_lb = max(certified_lb, m_final.leash)
        # Upper and lower bounds agree: mw is optimal.
        if (eq(m_final.leash, mw.leash, tolerance))
            return set_bounds!(mw, certified_lb)
        end

        # Not certified yet: smaller slack and a finer approximation.
        divisor = divisor * FC_SLACK_GROWTH
        aprx_refinement = 1.0 + ((aprx_refinement - 1.0) / FC_APPROX_SHRINK)
        # Close enough: the leash is within a factor 1 + tolerance of the lower bound.
        if f_accept_approx && ((1.0 + tolerance) * m_final.leash > mw.leash)
            return set_bounds!(mw, certified_lb)
        end
    end
end
