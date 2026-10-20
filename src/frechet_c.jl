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


"""
    frechet_c_approx

Approximates the continuous Frechet distance between the two input
curves. Returns a monotone morphing realizing it.

# Arguments

- `approx` : The output morhing has Frechet distance <= approx*optimal.

Importantly, approx can be larger than 2, if you want a rough
approximation.

The quality of approximation is available at ret.ratio. Thus,
ret.leash/ret.ratio is a lower bound on the Frechet distance, while
ret.leash is an upper bound.

"""
function frechet_c_approx(
    poly_a::PntSeq{N,T},
    poly_b::PntSeq{N,T},
    approx::Real,
)::Morphing{N,T} where {N,T}
    check_curves(poly_a, poly_b, "frechet_c_approx"; min_vertices = 2)
    approx = check_approx(approx, "frechet_c_approx")

    f_debug::Bool = false

    if f_debug
        println("\n\n\n============================XXX")
    end
    t_approx::Float64 = approx
    @assert((cardin(poly_a) > 1) && (cardin(poly_b) > 1))
    @assert(approx > 1.0)

    f_do_one_round::Bool = true
    d = frechet_dist_upper_bound(poly_a, poly_b)

    f_debug && println("DIST :", d)

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
        f_debug && println("AAAAAAAAAAAAAa")
        while ((r >= (d / (approx + 4.0))) || f_do_one_round)
            r = r / 2.0
            f_debug && println("BBBBBB r=", r)
            P, p_indices = Polygon_simplify_ext(poly_a, r)
            Q, q_indices = Polygon_simplify_ext(poly_b, r)

            if ((2 * cardin(P) > cardin(poly_a)) && (2 * cardin(Q) > cardin(poly_b)))
                f_debug && println("DDDD")
                f_debug && println("VVV Calling frechet_mono_via_refinement")
                m = frechet_mono_via_refinement(poly_a, poly_b)[1]
                m.ratio = 1.0
                m.lower_bound = m.leash
                #                exit( -1 ); # DDDD
                return m
            end

            @assert((cardin(P) > 1) && (cardin(Q) > 1))

            f_debug && println("cardins: ", cardin(P), " , ", cardin(Q))
            t_approx = 3.0/4.0 + t_approx / 4.0

            m = frechet_mono_via_refinement(P, Q, t_approx)[1]
            d = m.leash  #        frechet_ve_r_compute_dist( P, Q )
            f_debug && println("d : ", d)
            if (fp_equal(d, 0.0))
                f_debug && println("Breaking...")
                break
            end
            f_debug && println("ddd = ", d)
            f_debug && println("after")

            # d - 2r: Lower bound on the real Frechet distance.
            if f_debug
                println("d: ", d, "  r: ", r, "    #P: ", cardin(P), " Q: ", cardin(Q))
            end
            f_do_one_round = false
        end

        #mm = Morphing_monotonize( m );
        #        println( "mono a" );
        #m_p = frechet_ve_r_mono_compute( poly_a, P )
        #m_p = frechet_ve_r_mono_compute( poly_a, P )
        f_debug && println("=== f_c_dist( poly_a, P )")
        m_p = frechet_c_mono_approx_subcurve(poly_a, P, p_indices)[1]
        #        println( "mono b" );

        #m_q = frechet_ve_r_mono_compute( Q, poly_b )
        f_debug && println("f_c_dist( poly_b, Q )")
        m_q = frechet_c_mono_approx_subcurve(poly_b, Q, q_indices)[1]
        Morphing_swap_sides!(m_q)

        #err = max( m_p.leash, m_q.leash );
        #println( "Error: ", err );

        # m = m_p * mm * m_q
        f_debug && println("Combine m_p, m")
        mmu = Morphing_monotonize(Morphing_combine(m_p, m))
        f_debug && println("Combine mmu,m_q")
        mm_out = Morphing_combine(mmu, m_q)

        f_debug && println("Verification...")
        f_debug && Morphing_verify_valid(m_p)
        f_debug && Morphing_verify_valid(m)
        f_debug && Morphing_verify_valid(m_q)
        f_debug && Morphing_verify_valid(mmu)

        total_err = m_p.leash + m_q.leash
        lbx = (m.leash / t_approx) - total_err
        if f_debug
            println("total_error : ", total_err)
            println("lbx         : ", lbx)
        end
        if (total_err == 0.0) && (m.leash == 0.0)
            ratio = 1.0
        else
            if (lbx <= 0.0)
                ratio = 2.0 + approx # Very bad.
            else
                ratio = mm_out.leash / lbx
            end
        end
        f_debug && println("Ratio :", ratio)
        if ratio <= approx
            if f_debug
                println("-----------------------------------")
                println("Leash out: ", mm_out.leash)
                println("Leash lower bound: ", d - 2r)
                println("approx ratio          : ", ratio)
                println("Required approximation: ", 1.0 + approx)
                println(
                    "||P||: ",
                    fp_geom_length(poly_a),
                    "  ||Q||: ",
                    fp_geom_length(poly_b),
                )
            end
            mm_out.ratio = ratio
            mm_out.lower_bound = m.leash / ratio
            return mm_out

            #            break;
        end
        f_do_one_round = true
        if f_debug
            println("Another final round? ", ratio)
        end
    end
end


function count_below_zero(arr)
    count = 0
    for x in arr
        if x <= 0.0
            count = count + 1
        end
    end
    return count
end


function Polygon_simplify_radii_ext(P::PntSeq{N,T}, r::Vector{T}) where {N,T}
    f_exact::Bool = false
    #f_debug  &&  println( "#P :", cardin( P ) );
    #f_debug  &&  println( "#r :", length( r ) );

    PS, p_indices = Polygon_simplify_radii(P, r)
    #println( "@!@! |PS|: ", cardin( PS ) );
    @assert(length(p_indices) == cardin(PS))
    if (cardin(PS) <= (cardin(P) / 2))
        @assert(cardin(PS) == length(p_indices))
        return PS, p_indices, false
    end

    p_indices = Vector{Int64}()
    for i = 1:cardin(P)
        push!(p_indices, i)
    end
    #println( "______________________________" );
    @assert(cardin(P) == length(p_indices))
    return P, p_indices, true
end


function propogate_mins(qz, rounds::Int64 = 1)
    a = deepcopy(qz)
    for r = 1:rounds
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
    pl, ql = Morphing_extract_vertex_radii(m)

    lower_bound = m.lower_bound

    pz = ((lower_bound * ones(length(pl))) - pl) / factor
    qz = ((lower_bound * ones(length(ql))) - ql) / factor

    propogate_mins(pz, 1)
    propogate_mins(qz, 1)

    #println( pz );

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
frechet_continous_approx instead.

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

To really ensure converges, the monotone distance computed between the
simplification is computed using refinement, so tha the ve_r distance

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
    f_debug = false

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

    f_debug && println("P#", cardin(P))
    f_debug && println("Q#", cardin(Q))

    m_rough = frechet_c_approx(P, Q, aprx_refinement)
    if (m_rough.leash == 0.0)
        return m_rough
    end

    ratio_2 = m_rough.ratio
    if f_debug
        println("RATIO 2: ", m_rough.ratio, "\nleash: ", m_rough.leash)
    end

    len_a = fp_geom_length(P)
    len_b = fp_geom_length(Q)
    f_debug && println("#", cardin(P))

    approx = min(1.0 + (len_a + len_b) / (100.0 * m_rough.leash), 1.1)
    if f_debug
        println(
            "\n\nRatio         : ",
            ratio_2,
            "\n\nlen_a         : ",
            len_a,
            "\nlen_b         : ",
            len_b,
            "\n2approx_leash : ",
            m_rough.leash,
            "\n\n",
            approx,
            " approximation...",
        )
    end

    m = m_rough
    if m_rough.ratio > approx
        m = frechet_c_approx(P, Q, approx)
    end
    ratio = m.ratio

    if (f_debug)
        println("Approximation computed...")
        println("m.leash : ", m.leash)
    end

    # Okay - we computed a good enough approximation. 
    # We next extract the leases distances at the vertices...
    pl, ql = Morphing_extract_vertex_radii(m)

    lower_bound = m.leash / ratio
    if f_debug
        println("Before while loop...")
        println("LOWER BOUND: ", lower_bound)
        println("upper BOUND: ", m.leash, "\n")
    end

    ### Now we start looping refining if we need to improve th eapproxiamtion.
    round::Int64 = 0
    while true
        round = round + 1
        check_rounds(round, "frechet_c_compute")
        aprx_refinement = max(aprx_refinement, min_approx)

        f_debug && println(
            "---------------------------\nFactor: ",
            factor,
            "  Approx : ",
            aprx_refinement,
        )

        # Computes the slack at the vertices. The larger tha slack the more one can move a 
        # vertex without changing the Frechet distance. Negatie slack means not to move the 
        # vertices.
        pz = ((lower_bound * ones(length(pl))) - pl) / factor
        qz = ((lower_bound * ones(length(ql))) - ql) / factor

        # Propogate the do not move by doing a local min convolution. The current code is 
        # not efficient.
        propogate_mins(pz, 8)
        propogate_mins(qz, 8)

        if f_debug
            println("pz count <0: ", count_below_zero(pz), " / ", length(pz))
            println("qz count <0: ", count_below_zero(qz), " / ", length(qz))
        end
        f_PS_exact::Bool = false
        f_QS_exact::Bool = false

        # Simplify the polygons using the slack computed above per vertex.
        PS, p_indices, f_PS_exact = Polygon_simplify_radii_ext(P, pz)
        QS, q_indices, f_QS_exact = Polygon_simplify_radii_ext(Q, qz)

        if f_debug
            println("|PS| = ", cardin(PS), "  |P| = ", cardin(P))
            println("|QS| = ", cardin(QS), "  |Q| = ", cardin(Q))
        end

        # If the PS and QS are really large, this is a waste of
        # time. Just compute the Frechet distance using
        # monotonization.
        if ((2 * cardin(PS) > cardin(P)) && (2 * cardin(QS) > cardin(Q)))
            f_debug && println("Simpl. not profitable, doing mono")
            return frechet_mono_via_refinement(P, Q, 1.0 + tolerance)[1]
        end

        t_aprx = max(tolerance + 1.0, aprx_refinement)
        f_debug && println("f_m_via_refine( PS, QS )\nt_aprx : ", t_aprx)

        m_mid, _f_exact, PSR, QSR = frechet_mono_via_refinement(PS, QS, t_aprx)

        f_debug && println("frechet mono via refinement computed")

        # BUG FIX: Used the mono computation instead of the subcurve code,
        # which is much faster...
        #println( "_#PS: ", cardin( PS ) );
        m_a, offsets_a = frechet_c_mono_approx_subcurve(P, PS, p_indices)
        mmu = Morphing_combine(m_a, m_mid)
        m_b, offsets_b = frechet_c_mono_approx_subcurve(Q, QS, q_indices)

        Morphing_swap_sides!(m_b)

        mw = Morphing_combine(mmu, m_b)

        f_debug && println("eq?  ", m_mid.leash, " = ", mw.leash)

        # is there is NO hope we got the optimal solution?
        if (! eq(m_mid.leash, mw.leash, tolerance))
            lb = mw.leash - 2.0*m_a.leash - 2.0*m_b.leash
            lower_bound = max(lower_bound, lb)
            if f_debug
                println("!!! Ratio        : ", fp_ratio(m_mid.leash, mw.leash))
                println("!!! m_mid.leash  : ", m_mid.leash)
                println("!!! mw.leash     : ", mw.leash)
                println("!!! Diff         : ", m_mid.leash - mw.leash)
                println("!!! lower_bound  : ", lower_bound)
                println("!!! LB           : ", lb)
            end
            factor = factor * factor_scale
            aprx_refinement = 1.0 + (aprx_refinement - 1.0) / approx_scale
            f_debug &&
                println("1 aprx_refinement: ", aprx_refinement, "  scale :", approx_scale)

            pl, ql = Morphing_extract_vertex_radii(mw)
            continue
        end

        f_debug && println("Extracting offsets...")

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
        f_debug && println("*** lower_bound.leash  : ", m_final.leash)
        f_debug && println("*** mw     .leash      : ", mw.leash)

        if (eq(m_final.leash, mw.leash, tolerance))
            f_debug && println("Return from frechet_c_compute")
            return mw
        end

        # WE have to try harder - making slack at vertices smaller, and approximation smaller...
        factor = factor * factor_scale
        aprx_refinement = 1.0 + ((aprx_refinement - 1.0) / approx_scale)
        if f_debug
            println("m_a.leash      : ", m_a.leash)
            println("m_b.leash      : ", m_b.leash)
            println("m_mid.leash    : ", m_mid.leash)
            println("MW Leash length: ", mw.leash)
            println("m_final.leash  : ", m_final.leash)
        end
        if f_accept_approx && (1.00001 * m_final.leash > mw.leash)
            f_debug && println("MW Leash length: ", mw.leash)
            f_debug && println("m_final.leash  : ", m_final.leash)
            f_debug && println("Return (2) from frechet_c_compute")
            return mw
        end
    end  # While loop end
end

