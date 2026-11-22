###########################################################################
# Compute the morphing src(v) -> trg(u): That is u(v(t))
###########################################################################

"""
    frechet_dist_upper_bound

Returns a rough upper bound on the Frechet distance between the two
curves. This upper bound is on the continuous distance. No guarenteee
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

Computing the ver frechet distance, and then monotonize it. No
guarentee as far as the quality of the distance this realizes...

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
    m.ratio = (m.lower_bound > 0.0) ? m.leash / m.lower_bound : (m.leash == 0.0 ? 1.0 : Inf)
    return m
end


# Upper/lower ratio certified for a morphing built from simplifications with
# total error total_err.
function approx_ratio_bound(m_leash, t_approx, total_err, out_leash, approx)
    lbx = (m_leash / t_approx) - total_err
    if (total_err == 0.0) && (m_leash == 0.0)
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

- `approx` : The output morhing has Frechet distance <= approx*optimal.

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
    for i in 1:cardin(P)
        push!(p_indices, i)
    end
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
distance senstiave simplification of the curves. It then compute the
monotone fr_ve_r distance between the simplified curves, and it
combine it to get a distance between the two original cuves. It makre
sure the answers are the same, otherwise, it repeates with a finer
simplification/approximation till they are equal.

Finally, the algorithm uses the fr_ve_r_with_offests distance between
the two simplified curves to comptue a lower bound, and make sure this
is equal to the Frechet distance computed. If they are equal, then the
upper/lower bounds on the Frechet distance of the two curves are the
same, which implies that the computed distance is indeed the desired
Frechet distance.

# More details

To ensure convergence, the monotone distance between the simplified curves is
computed by refinement (`frechet_mono_via_refinement`), so that it approaches their
VE-retractable distance.

# Arguments
- f_accept_approx: accept a good enough approximation as good enough...

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

    # The parameters that can be fine-tuned....
    aprx_refinement::Float64 = 2.0
    factor::Float64 = 4.0
    factor_scale::Float64 = 1.3
    approx_scale::Float64 = 10.0

    min_approx::Float64 = 1.00000000001


    m_rough = frechet_c_approx(P, Q, aprx_refinement)
    if (m_rough.leash == 0.0)
        return m_rough
    end


    len_a = geom_length(P)
    len_b = geom_length(Q)

    approx = min(1.0 + (len_a + len_b) / (100.0 * m_rough.leash), 1.1)
    # Curves far apart relative to their length round this to exactly 1.
    if approx <= 1.0
        approx = min_approx
    end

    m = m_rough
    if m_rough.ratio > approx
        m = frechet_c_approx(P, Q, approx)
    end
    ratio = m.ratio


    # Okay - we computed a good enough approximation. 
    # We next extract the leases distances at the vertices...
    pl, ql = Morphing_extract_vertex_radii(m)

    lower_bound = m.leash / ratio
    # A certified lower bound, reported in the result; lower_bound above only steers
    # the simplification.
    certified_lb = m.lower_bound

    ### Now we start looping refining if we need to improve th eapproxiamtion.
    round::Int64 = 0
    while true
        round = round + 1
        check_rounds(round, "frechet_c_compute")
        aprx_refinement = max(aprx_refinement, min_approx)


        # Computes the slack at the vertices. The larger tha slack the more one can move a 
        # vertex without changing the Frechet distance. Negatie slack means not to move the 
        # vertices.
        pz = ((lower_bound * ones(length(pl))) - pl) / factor
        qz = ((lower_bound * ones(length(ql))) - ql) / factor

        # Propogate the do not move by doing a local min convolution. The current code is 
        # not efficient.
        propagate_mins!(pz, 8)
        propagate_mins!(qz, 8)

        f_PS_exact::Bool = false
        f_QS_exact::Bool = false

        # Simplify the polygons using the slack computed above per vertex.
        PS, p_indices, f_PS_exact = Polygon_simplify_radii_ext(P, pz)
        QS, q_indices, f_QS_exact = Polygon_simplify_radii_ext(Q, qz)


        # If the PS and QS are really large, this is a waste of
        # time. Just compute the Frechet distance using
        # monotonization.
        if ((2 * cardin(PS) > cardin(P)) && (2 * cardin(QS) > cardin(Q)))
            m = frechet_mono_via_refinement(P, Q, 1.0 + tolerance)[1]
            return set_bounds!(m, max(certified_lb, m.lower_bound))
        end

        t_aprx = max(tolerance + 1.0, aprx_refinement)

        m_mid, _f_exact, PSR, QSR = frechet_mono_via_refinement(PS, QS, t_aprx)


        # BUG FIX: Used the mono computation instead of the subcurve code,
        # which is much faster...
        m_a, offsets_a = frechet_c_mono_approx_subcurve(P, PS, p_indices)
        mmu = Morphing_combine(m_a, m_mid)
        m_b, offsets_b = frechet_c_mono_approx_subcurve(Q, QS, q_indices)

        Morphing_swap_sides!(m_b)

        mw = Morphing_combine(mmu, m_b)
        # m_mid.lower_bound <= d(PS, QS) <= d(P, Q) + m_a.leash + m_b.leash.
        certified_lb = max(certified_lb, m_mid.lower_bound - m_a.leash - m_b.leash)


        # is there is NO hope we got the optimal solution?
        if (! eq(m_mid.leash, mw.leash, tolerance))
            lb = mw.leash - 2.0*m_a.leash - 2.0*m_b.leash
            lower_bound = max(lower_bound, lb)
            factor = factor * factor_scale
            aprx_refinement = 1.0 + (aprx_refinement - 1.0) / approx_scale

            pl, ql = Morphing_extract_vertex_radii(mw)
            continue
        end


        # Now we compute the distance, with offsets...
        #PS_offs = offsets_a;#Morphing_extract_offsets( m_a )[2]
        #QS_offs = offsets_b;#Morphing_extract_offsets( m_b )[1]

        map_PS = get_refinement_map(PSR, PS)
        map_QS = get_refinement_map(QSR, QS)

        PS_offs = offsets_copy_map(map_PS, offsets_a)
        QS_offs = offsets_copy_map(map_QS, offsets_b)
        @assert(length(PS_offs) == cardin(PSR))
        @assert(length(QS_offs) == cardin(QSR))

        # very important! We have to use the refined curves PSR/QSR to monotonicity 
        # copies...
        m_final = frechet_ve_r_compute_ext(PSR, QSR, PS_offs, QS_offs, true)

        # The VE-retractable distance with offsets is a lower bound on d(P, Q).
        certified_lb = max(certified_lb, m_final.leash)
        if (eq(m_final.leash, mw.leash, tolerance))
            return set_bounds!(mw, certified_lb)
        end

        # WE have to try harder - making slack at vertices smaller, and approximation smaller...
        factor = factor * factor_scale
        aprx_refinement = 1.0 + ((aprx_refinement - 1.0) / approx_scale)
        if f_accept_approx && (1.00001 * m_final.leash > mw.leash)
            return set_bounds!(mw, certified_lb)
        end
    end  # While loop end
end
