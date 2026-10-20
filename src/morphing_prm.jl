"""
    extract_param_inner

Takes the sequence of events (i.e., sequence of points along the
curve), and outputs a corresponding sequence of real numbers, where
the ith nubmer if the distance of the ith point in sequence from the
beginning of the curve. Here, distance means the total distance of the
subcurve from the start point, to the current point.
"""
function extract_param_inner(P::PntSeq{N,T}, V::Vector{EventPoint{N,T}}) where {N,T}
    lens::Vector{Float64} = fp_prefix_lengths(P)
    out = Vector{Float64}()
    check_times(V)

    n = cardin(P)
    for i in eachindex(V)
        ev::EventPoint{N,T} = V[i]
        if ev.type == PT_VERTEX
            push!(out, lens[ev.i])
            continue
        end
        #println( "ev.t: ", ev.t, "  i: ", i,"  len(V): ", length( V )  );
        curr = lens[ev.i]
        next = lens[min(n, ev.i + 1)]
        push!(out, curr + ev.t * (next - curr))
    end

    return out
end


"""
    Morphing_extract_prm

    A parameterization is a polygonal curve that starts at (0,0) and
    end at (m,n).  The polygonal curve either have positive slope edge,
    or vertical or horizontal edges. It can be thought of as a
    piecewise linear function from [0,m] to [0,n]. Here m and n are the
    lengths of the two given polygons of P and Q, respectively.

"""
function Morphing_extract_prm(m::Morphing{N,T})::PntSeq2F where {N,T}
    P::PntSeq{N,T} = m.P
    Q::PntSeq{N,T} = m.Q
    peout::Vector{EventPoint{N,T}} = m.pes
    qeout::Vector{EventPoint{N,T}} = m.qes

    check_times(m.pes)

    check_times(m.qes)

    # pps and qps are two
    pps::Vector{Float64} = extract_param_inner(P, peout)
    qps::Vector{Float64} = extract_param_inner(Q, qeout)

    @assert(length(pps) == length(qps))

    out = PntSeq2F()
    for i in eachindex(pps)
        push!(out, point(pps[i], qps[i]))
    end
    return out
end


function eval_pl_func_on_dim(
    p::Point{N,T},
    q::Point{N,T},
    val::Float64,
    d::Int64,
) where {N,T}
    @assert(p[d] <= q[d])
    t = (val - p[d]) / (q[d] - p[d])

    return fp_convex_comb(p, q, t)
end

function eval_pl_func(p::Point2F, q::Point2F, val::Float64)
    o = eval_pl_func_on_dim(p, q, val, 1)
    return o[2]
end
function eval_inv_pl_func(p::Point2F, q::Point2F, val::Float64)
    o = eval_pl_func_on_dim(p, q, val, 2)
    return o[1]
end

function check_monotone_top(out::PntSeq2F)
    l = cardin(out)
    if l < 2
        return
    end
    p = out[l-1]
    q = out[l]
    factor::Float64 = 1.0001
    if (p[1] > factor*q[1]) || (p[2] > factor*q[2])
        throw(
            ArgumentError(
                "parameterizations are not monotone ($(p[1]) > $(q[1]) or $(p[2]) > $(q[2])); " *
                "Morphing_combine needs monotone morphings",
            ),
        )
    end
end



##############################################################3
# parameterization_combine
# Given two parameterization f, g, compute f(g( . ))
###############################################################
function parameterization_combine(f::PntSeq2F, g::PntSeq2F)::PntSeq2F
    fp_equal(last(g)[2], last(f)[1]) || throw(
        ArgumentError(
            "parameterization_combine: the parameterizations do not compose " *
            "($(last(g)[2]) != $(last(f)[1]))",
        ),
    )
    idf::Int64 = 1
    idg::Int64 = 1
    out::PntSeq2F = PntSeq2F()

    l_f = cardin(f)
    l_g = cardin(g)

    @assert(l_f > 1)
    @assert(l_g > 1)


    # x = dim(g,1)
    # y = dim(g,2) = dim(f, 1)
    # z = dim(f,2)
    while (true)
        check_monotone_top(out)

        # The two points under consideration, on the y axis
        yf = f[idf][1]
        yg = g[idg][2]
        f_equal_yf_yg::Bool = fp_equal(yf, yg)
        if (f_equal_yf_yg) && (idf == l_f) && (idg == l_g)
            push!(out, point(g[idg][1], f[idf][2]))
            break
        end
        if (f_equal_yf_yg && (idf < l_f) && (f[idf+1][1] == yf))
            push!(out, point(g[idg][1], f[idf][2]))
            idf = idf + 1
            continue
        end
        if (f_equal_yf_yg && (idg < l_g) && (fp_equal(g[idg+1][2], yg)))
            push!(out, point(g[idg][1], f[idf][2]))
            idg = idg + 1
            continue
        end
        if f_equal_yf_yg
            push!(out, point(g[idg][1], f[idf][2]))
            idf = min(idf + 1, l_f)
            idg = min(idg + 1, l_g)
            continue
        end
        if (yf < yg)
            #            println( "idf :", idf  );
            #            println( "idg :", idg );
            #            println( "yf :" ,yf );
            #            println( "yg :" ,yg );
            #            println( g );

            # compute g( yf )...
            xg = eval_inv_pl_func(g[idg-1], g[idg], yf)
            #        println( "xg: ", xg );

            # A bit of a hack... Because of floating point errors,
            # things can go a bit backward, which we really should not
            # allow.
            if (xg < g[idg-1][1])
                xg = g[idg-1][1]
            end
            push!(out, point(xg, f[idf][2]))
            idf = min(idf + 1, l_f)
            continue
        end
        if (yf > yg)
            zf = eval_pl_func(f[idf-1], f[idf], yg)
            # A bit of a hack again...
            if (zf > f[idf][2])
                zf = f[idf][2]
            end
            if (zf < f[idf-1][2])
                zf = f[idf-1][2]
            end
            push!(out, point(g[idg][1], zf))
            idg = min(idg + 1, l_g)
            continue
        end
        @assert(false)
    end

    # We monotononize the output because of floating point error...
    # Should really just get rid of tiny flatting point jitters.
    curr::Point2F = deepcopy(out[1])
    for i = 2:cardin(out)
        curr = point_max(out[i], curr)
        out[i] = curr
        #=for  j in 1:2
            curr[ j ] = max( curr[ j ], out[ i ][ j ] )
            out[ i ][ j ] = curr[ j ];
        end=#
    end

    return out
end

"""
    get_point

Returns the point along P that is in distance pos from the begining of
P ( lenfgth. This is the inner function, where we also provide the
edge it lies on, and the prefix sums (lens) of the lengths of the
edgss of P.

"""
function get_point(
    P::PntSeq{N,T},
    lens::Vector{Float64},
    i::Int64,
    pos::Float64,
) where {N,T}

    if (i == length(lens))
        return last(P), 1
    end
    edge_len = lens[i+1] - lens[i]
    t = (pos - lens[i]) / edge_len

    if (0 > t > -0.00000001)
        t = 0.0
    end
    if (1.0 < t < 1.000000001)
        t = 1.0
    end

    #    println( "inside... 3" );
    if (t == 0)
        return P[i], 0.00
    end
    if (t == 1)
        return P[i+1], 1.0
    end
    p = fp_convex_comb(P[i], P[i+1], t)
    #    println( "inside... 4" );
    return p, t
end

function event_sequences_extract(prm::PntSeq2F, P::PntSeq{N,T}, Q::PntSeq{N,T}) where {N,T}

    i_p::Int64 = 1
    i_q::Int64 = 1
    pes = Vector{EventPoint{N,T}}() # P event sequence
    qes = Vector{EventPoint{N,T}}() # Q event sequence

    len_q = cardin(Q)
    len_p = cardin(P)

    lp::Vector{Float64} = fp_prefix_lengths(P)
    lq::Vector{Float64} = fp_prefix_lengths(Q)

    for i = 1:(cardin(prm)-1)
        #        println( "before?? A" );
        curr = prm[i]
        p_loc = curr[1]

        while (i_p < length(lp)) && (p_loc >= lp[i_p+1])
            i_p = i_p + 1
        end
        (lp[i_p] <= p_loc) || error(
            "event_sequences_extract: parameter $p_loc precedes the vertex position " *
            "$(lp[i_p]) (vertex $i_p of $(length(lp)))",
        )
        while (i_q < length(lq)) && (curr[2] >= lq[i_q+1])
            i_q = i_q + 1
        end

        #t_p::Float64;
        pcurr, t_p = get_point(P, lp, i_p, p_loc)
        @assert(0 <= t_p && t_p <= 1.0)
        #@assert( ! isNaN(pcurr) );
        #        println( "get_point Q ..." );
        #get_point( Q, lq, i_q, curr[ 2 ] );
        #        println( "get_point Q!!!!!!!!!! ..." );
        qcurr, t_q = get_point(Q, lq, i_q, curr[2])
        #@assert( ! isNaN( qcurr ) );

        (0.0 <= t_q <= 1.0) ||
            error("event_sequences_extract: the parameter t_q = $t_q is outside [0, 1]")
        #        println( "here? " );

        if (t_p == 0.0) || ((t_p == 1.0) && (i_p == len_p))
            push!(pes, EventPoint(deepcopy(pcurr), i_p, PT_VERTEX, 0.0, 0.0))
        else
            push!(pes, EventPoint(deepcopy(pcurr), i_p, PT_ON_EDGE, T(t_p), 0.0))
        end

        #println( pcurr );
        #println( qcurr );
        if (t_q == 0.0) || ((t_q == 1) && (i_q == len_q))
            push!(qes, EventPoint(deepcopy(qcurr), i_q, PT_VERTEX, 0.0, 0.0))
        else
            push!(qes, EventPoint(deepcopy(qcurr), i_q, PT_ON_EDGE, T(t_q), 0.0))
        end
        #       println( "here? B" );
    end

    @assert(!any(isnan, last(P)))
    @assert(!any(isnan, last(Q)))

    push!(pes, EventPoint(deepcopy(last(P)), cardin(P), PT_VERTEX, 0.0, 0.0))
    push!(qes, EventPoint(deepcopy(last(Q)), cardin(Q), PT_VERTEX, 0.0, 0.0))

    #=
    println( "@@@@@@@@@@@@@@@@@@@@222" );
    println( pes );
    println( "@@@@@@@@@@@@@@@@@@@@222" );
    println( qes );
    =#
    #=
    for  i in 1:length( qes )
        if ( ( pes[ i ].type == PT_ON_EDGE ) &&
            ( qes[ i ].type != PT_ON_EDGE ) )
            println( "double edge event?" );
            println( "i:", i );
            println( "qes[ i ].t:", qes[ i ].t );
            println( "pes[ i ].t:", pes[ i ].t );
            println( length( qes ) );
        end
        @assert( ( pes[ i ].type != PT_ON_EDGE )
                 ||  ( qes[ i ].type != PT_ON_EDGE ) );
    end
    =#

    return pes, qes
end

########################################################################
"""
    Morphing_from_prm( prm, P, Q )

    Construct a morphing from a parameterization of the two polygons P
    and Q.

"""
function Morphing_from_prm(prm::PntSeq2F, P::PntSeq{N,T}, Q::PntSeq{N,T}) where {N,T}
    pes, qes = event_sequences_extract(prm, P, Q)
    return Morphing_init(P, Q, pes, qes)
end


function compute_offsets(
    P::PntSeq{N,T},
    pes::Vector{EventPoint{N,T}},
    qes::Vector{EventPoint{N,T}},
) where {N,T}
    offs = zeros(Float64, cardin(P))

    for i in eachindex(pes)
        ev = pes[i]
        evx = qes[i]
        r = fp_dist(ev.p, evx.p)
        offs[ev.i] = max(offs[ev.i], r)
    end

    return offs
end

"""
    Morphing_extract_offsets
"""
function Morphing_extract_offsets(m::Morphing)
    p_offs = Vector{Float64}()
    q_offs = Vector{Float64}()


    p_offs = compute_offsets(m.P, m.pes, m.qes)
    q_offs = compute_offsets(m.Q, m.qes, m.pes)
    # P event sequence

    return p_offs, q_offs
end





