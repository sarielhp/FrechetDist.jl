#####################################################################
"""
    discrete_frechet_sample

Computes the discrete Frechet distance between the two curves.by
sampling them

It first sapmles the two curves rougly uniformly. n is supposed to
be the nubmer of vertices, by the optput might be a sample that is
slightly bigger, as the code tries to ensure the vertices are being
picked.

# Arguments

- `n`: number of vertices to add when "refining" the two curves. The
  number of vertices computed might be larger (but hopefully not much
  larger).

- f_lopt = true: Use standard discrete Frechet or the retractable
  version.

"""
#####################################################################
function discrete_frechet_sample(
    polya::PntSeq{N,T},
    polyb::PntSeq{N,T},
    n::Int64,
    f_lopt::Bool = true,
) where {N,T}
    #    polya,polyb = example_4()
    lena = fp_geom_length(polya)
    lenb = fp_geom_length(polyb)

    #    n::Int64 = 100;
    delta = (lena + lenb) / n

    na::Int64 = ceil(lena / delta) + 2
    nb::Int64 = ceil(lenb / delta) + 2

    qa = Polygon_sample_uniformly(polya, na)
    qb = Polygon_sample_uniformly(polyb, nb)

    if (f_lopt)
        morph = frechet_d_r_compute(qa, qb)
    else
        morph = frechet_d_compute(qa, qb)
    end
    #    output_morphing( morph, filename, polyb );

    return morph
end


"""
    frechet_c_mono_approx_subcurve

Approximates the Frechet distance between a curve (P) and subcurbe (Q). Here,
Q vertices are the vertices of P specifieid by p_indices. That is
Q[ i ] =P[ p_indices[ i ]].

"""
function frechet_c_mono_approx_subcurve(
    P::PntSeq{N,T},
    Q::PntSeq{N,T},
    p_indices::Vector{Int64},
) where {N,T}

    cardin(Q) == length(p_indices) || throw(
        ArgumentError(
            "frechet_c_mono_approx_subcurve: Q has $(cardin(Q)) vertices but " *
            "$(length(p_indices)) indices were given",
        ),
    )
    check_indices(p_indices, cardin(P), "frechet_c_mono_approx_subcurve")
    offsets = zeros(T, length(p_indices))

    pes = Vector{EventPoint{N,T}}()
    qes = Vector{EventPoint{N,T}}()

    for i in eachindexButLast(p_indices)
        pind = p_indices[i]
        pind_next = p_indices[i+1]
        push!(pes, EventPoint(P[pind], pind, PT_VERTEX, 0.0, 0.0))
        push!(qes, EventPoint(Q[i], i, PT_VERTEX, 0.0, 0.0))

        max_d::Float64 = 0
        max_t::Float64 = 0
        seg = Segment(Q[i], Q[i+1])
        for j = (pind+1):(pind_next-1)
            p = P[j]
            q = fp_nn_point(seg, p)
            new_t = fp_convex_coef(seg, q)
            if (new_t >= max_t)
                max_t = new_t
            else # new_t < max_t
                q = fp_at(seg, max_t)
            end
            max_d = max(max_d, fp_dist(P[j], q))
            push!(pes, EventPoint(P[j], j, PT_VERTEX, 0.0, 0.0))
            push!(qes, EventPoint(q, i, PT_ON_EDGE, max_t, 0.0))
        end
        offsets[i] = max(offsets[i], max_d)
        offsets[i+1] = max(offsets[i+1], max_d)
    end

    push!(pes, EventPoint(last(P), cardin(P), PT_VERTEX, 0.0, 0.0))
    push!(qes, EventPoint(last(Q), cardin(Q), PT_VERTEX, 0.0, 0.0))

    return Morphing_init(P, Q, pes, qes), offsets
end


"""
    add_points_along_seg

# Arguments

seg: t[1], t[2], .... , t[end]

     t[...]: Time stamps for points on the segment (all between 0 and 1)

# Details

Add points along seg at time stamps specified by t, and add points in
the middle between them (i.e., refining segment). These new points are
added the output polygon pout.

Note, that one can set n_mid_points to be larger than 1, but it does
not seem to be necessary/benefitial from the experimentation we did.

Output:
     Sorted( t[1]/2, t[1], (t[1]+t[2])/2, t[2],...)
"""
function add_points_along_seg(
    pout::PntSeq{N,T},
    seg::Segment{N,T},
    times::Vector{Float64},
    n_mid_points::Int64 = 1,
) where {N,T}
    out_times = Vector{Float64}()

    function push_time(tm::Float64)
        if (tm <= 0.0) || (tm >= 1.0)
            return
        end

        push!(out_times, tm)
    end

    sort!(times)
    unique!(times)
    times_length = length(times)

    if isempty(times)
        return
    end

    if (n_mid_points > 0)
        push_time(times[1] / 2.0)
    end

    for i = 1:length(times)
        t = times[i]
        push_time(t)
        if (i < times_length) && (n_mid_points > 0)
            K = n_mid_points + 1
            #println( "K : ", K );
            #if  ( times[ i + 1 ] > times[ i ] )  continue  end;
            for pos = 1:(K-1)
                coef::Float64 = pos / K
                tm = times[i] * (1.0 - coef) + times[i+1] * coef
                push_time(tm)
            end
        end
    end
    if (n_mid_points > 0)
        push_time((1.0 + last(times)) / 2.0)
    end

    if isempty(out_times)
        return
    end

    sort!(out_times)
    unique!(out_times)

    prev::Float64 = -1.0
    dlt = 0.09
    for tm in out_times
        if ((tm <= 0.0) || (tm >= 1.0))
            continue
        end
        if (tm < (prev + dlt))
            continue
        end
        prev = tm
        push_distinct!(pout, fp_at(seg, tm))
    end
end


function is_monotone_inc(nums::Vector{Float64})
    len = length(nums)
    if (len <= 1)
        return true
    end
    for i = 1:(len-1)
        if nums[i] > nums[i+1]
            return false
        end
    end

    return true
end


################################################################
"""
    extract_refined_polygon( P, pes, num_points_to_add, delta )

    Refine edges along which pes (the sequence of matching points
    along P is not monotone, so that the returned polygon has
    vertices, which hopefully force the Frechet morphing computed in
    the next iteration to be monotone.

# Arguments

  `P`  : the polygo
  `pes`: Sequence of events long the polygon P
  `num_points_to_add`: how many points to add between points of
                       non-monotonicity. The default is 1. Having a
                       larger value speeds up the process of getting
                       rid of non-monotonicity, but the polygons grow
                       faster.
"""
function extract_refined_polygon(
    poly::PntSeq{N,T},
    s::Vector{EventPoint{N,T}},
    points_to_add::Int64 = 1,
    delta::Float64 = 0.0,
) where {N,T}
    pout = PntSeq{N,T}()
    opnts = Vector{EventPoint{N,T}}()

    len = length(s)
    push!(pout, s[1].p)

    i = 2
    while (i <= len)
        ep = s[i]
        if ep.type == PT_VERTEX
            push_distinct!(pout, s[i].p)
            i = i + 1
            continue
        end

        j = i
        loc = s[i].i
        while ((j < len) && (s[j+1].type == PT_ON_EDGE) && (s[j+1].i == loc))
            j = j + 1
        end

        # s[i:j]: The sequence of points on the same edge
        times = Vector{Float64}()
        push!(times, s[i].t)
        for k = (i+1):j
            if s[k].t != last(times)
                push!(times, s[k].t)
            end
        end

        lnx = fp_dist(poly[loc], poly[loc+1])

        if (is_monotone_inc(times) || (lnx < (delta / 2.0)))
            for k = i:j
                push_distinct!(pout, s[k].p)
            end
            i = j + 1
            continue
        end

        #=
        diff = 0.0;
        min_t = max_t = times[ 1 ];

        for  i in 2:length(times)
            max_t = max( times[ i ], max_t );
            if  ( times[ i ] < max_t )  &&  ( ( max_t - times[ i ] ) > diff )
                min_t = times[ i ];
                diff = max_t - times[ i ];
            end
        end
        =#
        #println( "diff: ", diff );
        #println( "min_t: ", min_t );

        #println( times );
        #exit( -1 );
        n_times::Int64 = points_to_add
        if (delta > 0.0)
            n_times = round(Int64, lnx / delta)
            n_times = max(min(n_times, 3), 1)
        end
        seg = Segment(poly[loc], poly[loc+1])
        add_points_along_seg(pout, seg, times, n_times)

        i = j + 1
    end

    return pout
end


function frechet_refinement(P::PntSeq{N,T}, Q::PntSeq{N,T}, m::Morphing{N,T}) where {N,T}
    poly_a_2 = extract_refined_polygon(P, m.pes, 1)
    poly_b_2 = extract_refined_polygon(Q, m.qes, 1)

    return poly_a_2, poly_b_2
end


function frechet_mono_via_refinement_ext(
    Pa::PntSeq{N,T},
    Qa::PntSeq{N,T},
    out::Vector{Morphing{N,T}},
    f_snapshots::Bool,
    approx::Real = 1.00001,
) where {N,T}
    check_curves(Pa, Qa, "frechet_mono_via_refinement_ext"; min_vertices = 2)
    approx = check_approx(approx, "frechet_mono_via_refinement_ext"; exceed_one = false)
    if (approx <= 1.0)
        approx = 1.0000001
        @assert(approx > 1.0)
    end
    f_debug::Bool = false

    mm::Morphing{N,T} = Morphing_empty(Pa, Qa)
    m::Morphing{N,T} = Morphing_empty(Pa, Qa)
    f_exact::Bool = false
    f_refinement::Bool = false
    P = Pa
    Q = Qa
    f_debug && println("fmvr before while...")
    f_debug && println("APPROX : ", approx)
    rounds::Int = 0
    while true
        rounds += 1
        check_rounds(rounds, "frechet_mono_via_refinement_ext")
        f_debug && println("\n<<<<<<< ")
        f_debug && println("Before ve_r_compute( ", cardin(P), ", ", cardin(Q))
        m = frechet_ve_r_compute(P, Q)
        f_debug && println(">>> After... m.leash : ", m.leash)
        mm = Morphing_monotonize(m)

        if (f_snapshots)
            push!(out, m)
        end

        fr_retract = m.leash
        fr_r_mono = mm.leash

        if (f_debug)
            println("ve_r            : ", m.leash)
            println("ve_r.mono       : ", fr_r_mono)
            println("ve_r*approx     : ", approx * fr_retract)
        end

        if (fr_r_mono == fr_retract)
            f_exact = true
            break
        end
        if (fr_r_mono <= approx * fr_retract)
            f_exact = false
            break
        end
        f_refinement = true

        poly_a_2 = extract_refined_polygon(P, m.pes, 1)
        poly_b_2 = extract_refined_polygon(Q, m.qes, 1)

        P = poly_a_2
        Q = poly_b_2
    end

    # We potentially need to flatten the morphing back to the original
    # curves...
    if ((cardin(mm.P) != cardin(Pa)) || (cardin(mm.Q) != cardin(Qa)))
        prm = Morphing_extract_prm(mm)
        mm = Morphing_from_prm(prm, Pa, Qa)
    end

    mm.lower_bound = m.leash

    f_debug && println("RETURNING!", @__FILE__, ":", @__LINE__)
    return mm, f_exact, P, Q
end


function frechet_mono_via_refinement_delta(
    Pa::PntSeq{N,T},
    Qa::PntSeq{N,T},
    delta::Float64,
    f_fix_prm::Bool = false,
) where {N,T}
    check_curves(Pa, Qa, "frechet_mono_via_refinement_delta"; min_vertices = 2)
    check_nonnegative(delta, "frechet_mono_via_refinement_delta", "delta")
    f_debug::Bool = false

    P = Pa
    Q = Qa
    #println( "ROUNDS\n" );
    local m, mm
    count::Int64 = 0
    for count = 0:64
        #println( count, "  ", length( P ), "  ", length( Q ) );
        #( count > 0 )  &&  println( "FVER round: ", count );

        #@time m_old = frechet_ve_r_compute( P, Q )
        m = FEVER_compute(P, Q)
        #=
        println( "---" );

        if  ( ! fp_equal( m_old.leash, m.leash ) )
            @assert( m_old.leash == m.leash );
        end
        =#

        mm = Morphing_monotonize(m)

        res = 3
        if (delta > 0.0)
            gap = (mm.leash - m.leash) / delta
            if (gap > 0.0)
                res = 1 + (max(1, round(Int64, (log(gap)/log(2)))))^2
            end
        else
            gap = delta
        end

        ((mm.leash - m.leash) <= delta) && break
        #println( count, ": X GAP: ", gap, "   res: ", res );
        #println( count, ":   res: ", res );

        #times::Int64 = 4;# min( round(Int64,
        #   (mm.leash - m.leash) / delta ), 7 );;
        # res was 3
        poly_a_2 = extract_refined_polygon(P, m.pes, res, delta/11)
        poly_b_2 = extract_refined_polygon(Q, m.qes, res, delta/11)

        #        println( "P len: ", length( P ) );
        #println( "P_a_2 len: ", length( poly_a_2 ) );

        P = poly_a_2
        Q = poly_b_2
    end

    # We potentially need to flatten the morphing back to the original
    # curves...
    if f_fix_prm && ((cardin(mm.P) != cardin(Pa)) || (cardin(mm.Q) != cardin(Qa)))
        prm = Morphing_extract_prm(mm)
        mm = Morphing_from_prm(prm, Pa, Qa)
    end

    mm.lower_bound = m.leash

    return mm, P, Q
end



"""
    frechet_mono_via_refinement( P, Q, approx )

    Computes the "true" monotone Frechet distance between P and Q,
    using the ve_r algorithm. It does refinement, to add vertices if
    needed to get monotonicity. Note, that there is an eps parameter -
    for real inputs you can set it quite small. In any case, in the
    worst case it only approximates the Frechet (monotone)
    distance. It returns a morphing, and a boolean that is true if the
    result is the true Frechet distance.

    Observe, that this function might be too slow if the two curves
    are huge, since it does not simplify them before computing the
    distance.

# Returns

Returns a 4-tuple  ( m, f, P, Q ):
`m`: Morphing realizing result.
`f`: Is distance returned is the exact continuous Frechet distance.
`P`: Refined first input polygon
`Q`: Refined second input polygon

'out' if specified returns the middle morphings used in computing the
      solution.
"""
function frechet_mono_via_refinement(
    Pa::PntSeq{N,T},
    Qa::PntSeq{N,T},
    approx::Real = 1.00001,
) where {N,T}
    check_curves(Pa, Qa, "frechet_mono_via_refinement"; min_vertices = 2)
    approx = check_approx(approx, "frechet_mono_via_refinement"; exceed_one = false)
    ot::Vector{Morphing{N,T}} = Vector{Morphing{N,T}}()
    return frechet_mono_via_refinement_ext(Pa, Qa, ot, false, approx)
end

