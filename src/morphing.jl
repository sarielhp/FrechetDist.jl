# Originally contributed by S. Har-Peled
# under MIT License

#---------------------------------------------------------
# Morphing
#
# A morphing is the matching between two polygons as computed by by
# the Frechet distance.
#---------------------------------------------------------

using Parameters
using DataStructures
using Printf



@enum FPointType begin
    PT_VERTEX = 1
    PT_ON_EDGE = 2
end

############################################################
# A vertex edge event descriptor.
#
# p: Location of the point being matched. Not necessarily a vetex of
#    the polygon.
# i: Vertex number in polygon/or edge number where p lies.
# type:
# t: Convex combination coefficient if p is on the edge. t=0 means its
#    the ith vertex, t=1 means it is on the i+1 vertex.
############################################################
#@with_kw
mutable struct EventPoint{N,T}
    p::Point{N,T}
    i::Int64;
    type::FPointType;
    t::Float64;
end

# leash = maximum length of an edge in the encoded matching.
#         (That is, the Frechet distance.)
"""
    Morphing

Encoding of a moprhing (i.e., matching) between two polygonal cuves.

"""
@with_kw mutable struct  Morphing{N,T}
    P::Polygon{N,T}
    Q::Polygon{N,T}
    pes::Vector{EventPoint{N,T}}; # P event sequence
    qes::Vector{EventPoint{N,T}}; # Q event sequence
    leash::Float64
#    leash_offsets::Float64
    iters::Int64
    ratio::Float64
    monotone_err::Float64
    sol_value::Float64
end

function   Morphing_init( P::Polygon{N,T}, Q::Polygon{N,T},
    pes::Vector{EventPoint{N,T}},
    qes::Vector{EventPoint{N,T}} ) where {N,T}

    @assert( length( pes ) == length( qes ) );
    m = Morphing( P, Q, pes, qes,
                  0.0, 0, 0.0, 0.0, 0.0 );
    Morphing_recompute_leash( m );

    return  m;
end


####################################################
# Compute for every vertex the maximum leash length used with this
# vertex.
####################################################

function   extract_vertex_leash( P::Polygon{N,T},
    pes::Vector{EventPoint{N,T}},
    qes::Vector{EventPoint{N,T}}
) where {N,T}
    P_len = cardin( P );
    len = length( pes );
    vl = zeros( Float64, P_len )
    for  i in  1:len
        ev = pes[ i ];
        evx = qes[ i ];
        if  ( ev.type == PT_VERTEX )
            vl[ ev.i ] = max( vl[ ev.i ], Dist( ev.p, evx.p ) );
        end
    end

    return  vl;
end


"""
    Morphing_extract_vertex_radii

    Computes for each polygon vertex, the length of the longest edge
    in the matching attached ot it. It is a cheap upper bound on the
    local Frechet distance for each vertex (and implicitly the
    attached edge).
"""
function   Morphing_extract_vertex_radii( m::Morphing{N,T} ) where {N,T}
    pl = extract_vertex_leash( m.P, m.pes, m.qes );
    ql = extract_vertex_leash( m.Q, m.qes, m.pes );
    return  pl, ql
end


function   Morphing_swap_sides!( m::Morphing{N,T} ) where {N,T}
    m.P, m.Q = m.Q, m.P
#    swap!( m.P, m.Q );
    #swap!( m.pes, m.qes );
    m.pes, m.qes = m.qes, m.pes;
end

"""
    Morphing_empty

Returns an empty morhping (i.e., a constructor).
"""
function   Morphing_empty( P::Polygon{N,T}, Q::Polygon{N,T} )  where {N,T}
    pes = Vector{EventPoint{N,T}}();
    qes = Vector{EventPoint{N,T}}()
    r::Float64 = -1;
    return  Morphing( P, Q, pes, qes, r, 0, 0.0, 0.0, 0.0 );
end


###################################################################
# Check if the matching is monotone. That is the endpoints of the
# matching on each curve are sorted, in the same order as the
# matching.
##################################################################
function  events_seq_is_monotone( s::Vector{EventPoint{N,T}}
) where {N,T}

    len = length( s );
    i::Int64 = 1;
    while  ( i <= len )
        ep = s[ i ];
        if  ep.type == PT_VERTEX
            i = i + 1;
            continue;
        end

        j = i
        loc = s[ i ].i;
        while  ( ( j < len )
            &&  ( s[ j + 1 ].type == PT_ON_EDGE )
            &&  ( s[ j + 1 ].i == loc ) )
            j = j + 1
        end

        # i:j is the sequence of edge-vertex events
        for  k  in i:(j-1)
            if  s[ k + 1 ].t < s[ k ].t
	        return  false;
            end
        end
        i = j + 1;
    end

    return  true;
end


##############################################################################
# Make the matching (which might be potentially not monotone) into a monotone
# one.
##############################################################################
function  events_seq_make_monotone( P::Polygon{N,T},
                                    s::Vector{EventPoint{N,T}} ) where {N,T}
    ns = Vector{EventPoint{N,T}}();

    len = length( s );
    i::Int64 = 1;
    delta = 0.0;
    while  ( i <= len )
        ep = s[ i ];
        if  ep.type == PT_VERTEX
            push!( ns, deepcopy( ep ) )
            i = i + 1;
            continue;
        end

        j = i
        loc = s[ i ].i;
        t = s[ i ].t
        while  ( ( j < len )
                 &&  ( s[ j + 1 ].type == PT_ON_EDGE )
                 &&  ( s[ j + 1 ].i == loc ) )
            nep = deepcopy( s[ j ] );
            t = max( t, nep.t );
            if  ( t > nep.t )
                ell = Dist( P[ nep.i ], P[ nep.i + 1 ] );
                delta = max( delta, ( t - nep.t ) * ell );
                nep.t = t;
                nep.p = convex_comb( P[ nep.i ], P[ nep.i + 1 ],
                                     t  );
            end
            push!( ns, nep );
            j = j + 1
        end

        nep = deepcopy( s[ j ] );
        if  ( t > nep.t )
            ell = Dist( P[ nep.i ], P[ nep.i + 1 ] );
            delta = max( delta, ( t - nep.t ) * ell );
            nep.t = t;
            nep.p = convex_comb( P[ nep.i ], P[ nep.i + 1 ],
                                 t );
        end
        push!( ns, nep );

        i = j + 1;
    end

    @assert( length( ns ) == length( s ) );
    return  ns, delta;
end

function  Morphing_recompute_leash( m::Morphing{N,T} ) where  {N,T}
    @assert( length( m.pes ) == length( m.qes ) );
    r = 0;
    for  i  in eachindex( m.pes )
        ell = Dist( m.pes[ i ].p, m.qes[ i ].p )
        if   ell > r
            r = ell;
        end
    end
    if  ( m.leash != r )   &&   ( m.leash != 0 )
        println( "old leash:" , m.leash );
        println( "new leash:" , r );
        exit( -1 );
    end
    m.leash = r;
#    println(  " m.leash: ", m.leash );
end

"""
    Morphing_as_polygons

Turns the morphing matching into a "true" matching, by creating two
polygons that their edges are directly matched.  The output polygons P
and Q will definitely have reapeated points.

"""
function  Morphing_as_polygons( m::Morphing{N,T} ) where  {N,T}
    P = Polygon{N,T}();
    Q = Polygon{N,T}();

    for  i  in eachindex( m.pes )
        push!( P, m.pes[ i ].p );
        push!( Q, m.qes[ i ].p );
    end

    return  P, Q
end


function   Morphing_is_monotone( m::Morphing{N,T} ) where {N,T}
    return ( events_seq_is_monotone( m.pes )
             && events_seq_is_monotone( m.qes ) )
end

##########################################################3
# Turn the morphing into monotone morphing...
##########################################################3
"""
    Morphing_monotonize

Turns a morphing into a monotone morphing, by simply not going back,
staying in place if necessary.
"""
function  Morphing_monotonize( m::Morphing{N,T} ) where {N,T}
    if  Morphing_is_monotone( m )
        return  deepcopy( m );
    end
    P = m.P;
    Q = m.Q;

    pes_new, da = events_seq_make_monotone( P, m.pes );
    qes_new, db = events_seq_make_monotone( Q, m.qes );


    m_out = Morphing_init( P, Q, pes_new, qes_new );

    m_out.monotone_err = max( da, db );

    return  m_out;
end



##########################################################################


function  check_times( V::Vector{EventPoint{N,T}} ) where {N,T}
    for  ev in V
        if  ( ev.t < 0.0 )  ||  (ev.t > 1.0 )
            println( "Event time is wrong? ", ev.t );
            exit( -1 );
        end
    end
end


"""
    Morphing_verify_valid

Does some minimal checks that the morphing is valid. Speciifcally,
check the times stemps of the events are valid.
"""
function  Morphing_verify_valid( m::Morphing{N,T} ) where {N,T}
    check_times( m.pes );
    check_times( m.qes );
end


"""
    extract_param_inner

Takes the sequence of events (i.e., sequence of points along the
curve), and outputs a corresponding sequence of real numbers, where
the ith nubmer if the distance of the ith point in sequence from the
beginning of the curve. Here, distance means the total distance of the
subcurve from the start point, to the current point.
"""
function  extract_param_inner( P::Polygon{N,T},
                               V::Vector{EventPoint{N,T}} ) where {N,T}
    lens::Vector{Float64} = Polygon_prefix_lengths( P )
    out = Vector{Float64}();
    check_times( V );

    n = cardin( P );
    for  i in eachindex( V )
        ev::EventPoint{N,T} = V[ i ];
        if  ev.type == PT_VERTEX
            push!( out, lens[ ev.i ] );
            continue;
        end
        #println( "ev.t: ", ev.t, "  i: ", i,"  len(V): ", length( V )  );
        curr = lens[ ev.i ];
        next = lens[ min( n, ev.i + 1 ) ]
        push!( out, curr + ev.t * (next - curr) )
    end

    return  out
end


"""
    Morphing_extract_prm

    A parameterization is a polygonal curve that starts at (0,0) and
    end at (m,n).  The polygonal curve ither have positive slope edge,
    or vertical or horizontla edges. It can be thought of as a
    piecewise linear function from [0,m] to [0,n]. Here m and n are the
    lengths of the two given polygons of P and Q, respectively.

"""
function  Morphing_extract_prm( m::Morphing{N,T} )::Polygon2F where {N,T}
    P::Polygon{N,T} = m.P;
    Q::Polygon{N,T} = m.Q;
    peout::Vector{EventPoint{N,T}} = m.pes;
    qeout::Vector{EventPoint{N,T}} = m.qes;

    check_times( m.pes );

    check_times( m.qes );

    # pps and qps are two
    pps::Vector{Float64} = extract_param_inner( P, peout );
    qps::Vector{Float64} = extract_param_inner( Q, qeout );

    @assert( length( pps ) == length( qps ) )

    out = Polygon{2,Float64}();
    for  i in eachindex(pps)
        push!( out, point( pps[ i ], qps[ i ] ) )
    end
    return  out;
end


function  eval_pl_func_on_dim( p::Point{N,T}, q::Point{N,T}, val::Float64,
                               d::Int64 )  where {N,T}
    @assert( p[ d ] <= q[ d ] );
    t = (val - p[d]) / (q[d] - p[d]);

    return  p * (1.0 - t) + q * t;
end

function  eval_pl_func( p::Point2F, q::Point2F, val::Float64 )
    o = eval_pl_func_on_dim( p, q, val, 1 );
    return  o[2];
end
function  eval_inv_pl_func( p::Point2F, q::Point2F, val::Float64 )
    o = eval_pl_func_on_dim( p, q, val, 2 );
    return  o[1];
end

function  check_monotone_top( out::Polygon2F )
    l = cardin( out )
    if  l < 2
        return;
    end
    p = out[ l - 1 ];
    q = out[ l ];
    factor::Float64 = 1.0001
    if  ( p[ 1 ] > factor*q[1] ) || ( p[2] > factor*q[2] )
        println( p[ 1 ], "  >  ", q[1] );
        println( p[ 2 ], "  >  ", q[2] );
        println( "Error (not monotone top): " );
        for i in 1:l
            println( i, " : ", out[ i ] );
        end
        #println( out );
        println( "Error (not monotone top): " );
        exit( -1 );
    end
end

function  floating_equal( a::Float64, b::Float64 )::Bool
    if   a == b
        return  true;
    end
    return  abs( a - b ) <= (0.0000001* (abs( a)  + abs(b) ))
end


function  floating_equal( a, b )::Bool
    return  a == b;
end

##############################################################3
# parameterization_combine
# Given two parameterization f, g, compute f(g( . ))
###############################################################
function   parameterization_combine( f::Polygon2F, g::Polygon2F )
    if  ( ! floating_equal( last( g )[2], last( f )[ 1 ] ) )
        println( last( g )[2], " != ",  last( f )[ 1 ] );
        @assert( floating_equal( last( g )[2], last( f )[ 1 ] ) )
    end
    idf::Int64 = 1;
    idg::Int64 = 1;
    out::Polygon2F = Polygon2F();

    l_f = cardin( f )
    l_g = cardin( g )

    @assert( l_f > 1 );
    @assert( l_g > 1 );


    # x = dim(g,1)
    # y = dim(g,2) = dim(f, 1)
    # z = dim(f,2)
    while (true)
        check_monotone_top( out );

        # The two points under consideration, on the y axis
        yf = f[ idf ][1];
        yg = g[ idg ][2];
        f_equal_yf_yg::Bool = floating_equal( yf, yg )
        if  ( f_equal_yf_yg )  &&  ( idf == l_f )  &&  ( idg == l_g )
            push!( out, point( g[ idg ][1], f[ idf ][2] ) )
            break;
        end
        if  ( f_equal_yf_yg   &&  ( idf < l_f )
              &&  ( f[ idf + 1 ][1] == yf ) )
            push!( out, point( g[ idg ][1], f[ idf ][2] ) )
            idf = idf + 1
            continue;
        end
        if  ( f_equal_yf_yg   &&  ( idg < l_g )
              &&  ( floating_equal( g[ idg + 1 ][2], yg ) ) )
            push!( out, point( g[ idg ][1], f[ idf ][2] ) )
            idg = idg + 1
            continue;
        end
        if  f_equal_yf_yg
            push!( out, point( g[ idg ][1], f[ idf ][2] ) )
            idf = min( idf + 1, l_f );
            idg = min( idg + 1, l_g );
            continue;
        end
        if  ( yf < yg )
#            println( "idf :", idf  );
#            println( "idg :", idg );
#            println( "yf :" ,yf );
#            println( "yg :" ,yg );
#            println( g );

            # compute g( yf )...
            xg = eval_inv_pl_func( g[idg - 1], g[ idg ], yf )
    #        println( "xg: ", xg );

            # A bit of a hack... Because of floating point errors,
            # things can go a bit backward, which we really should not
            # allow.
            if  ( xg < g[ idg - 1 ][ 1 ] )
                xg = g[ idg - 1 ][ 1 ]
            end
            push!( out, point( xg, f[ idf ][ 2 ] ) )
            idf = min( idf + 1, l_f );
            continue;
        end
        if  ( yf > yg )
            zf = eval_pl_func( f[idf - 1], f[ idf ], yg )
            # A bit of a hack again...
            if  ( zf > f[ idf ][ 2 ] )
                zf = f[ idf ][ 2 ]
            end
            if  ( zf < f[ idf - 1 ][ 2 ] )
                zf = f[ idf - 1 ][ 2 ]
            end
            push!( out, point( g[ idg ][ 1 ], zf ) )
            idg = min( idg + 1, l_g );
            continue;
        end
        @assert( false );
    end

    # We monotononize the output because of floating point error...
    # Should really just get rid of tiny flatting point jitters.
    curr::Point2F = deepcopy( out[ 1 ] );
    for  i in 2:cardin( out )
        for  j in 1:2
            curr[ j ] = max( curr[ j ], out[ i ][ j ] )
            out[ i ][ j ] = curr[ j ];
        end
    end

    return  out;
end

"""
    get_point

Returns the point along P that is in distance pos from the begining of
P ( lenfgth. This is the inner function, where we also provide the
edge it lies on, and the prefix sums (lens) of the lengths of the
edgss of P.

"""

function  get_point(
    P::Polygon{N,T}, lens::Vector{Float64},
    i::Int64, pos::Float64 ) where {N,T}

    if  ( i == length( lens ) )
        return  last( P ),1;
    end
    edge_len = lens[ i + 1 ] - lens[ i ];
    t = ( pos - lens[ i ] ) / edge_len;

    if  ( 0 > t > -0.00000001 )
        t = 0;
    end
    if  ( 1.0 < t < 1.000000001 )
        t = 1;
    end

#    println( "inside... 3" );
    if  ( t == 0 )
        return  P[ i ], .00;
    end
    if  ( t == 1 )
        return  P[ i + 1 ], 1.0;
    end
    p =  P[ i ] * (1.0 - t) + P[i + 1] * t;
#    println( "inside... 4" );
    return p, t;
end

function  event_sequences_extract( prm::Polygon2F, P::Polygon{N,T},
                                   Q::Polygon{N,T} ) where {N,T}

    i_p::Int64 = 1;
    i_q::Int64 = 1;
    pes = Vector{EventPoint{N,T}}(); # P event sequence
    qes = Vector{EventPoint{N,T}}(); # Q event sequence

    len_q = cardin( Q );
    len_p = cardin( P );

    lp::Vector{Float64} = Polygon_prefix_lengths( P )
    lq::Vector{Float64} = Polygon_prefix_lengths( Q )

    for  i in 1:( cardin( prm ) - 1 )
#        println( "before?? A" );
        curr::Point2F = prm[ i ];
        p_loc = curr[ 1 ];

        while  ( i_p < length( lp ) )  &&  ( p_loc >= lp[ i_p + 1 ] )
            i_p = i_p + 1;
        end
        if  ( ! ( lp[ i_p ] <= p_loc ) )
            println( "-----------------------" );
            println( "lp[ i_p - 1]  : ", lp[ i_p - 1 ] );
            println( "lp[ i_p ]     : ", lp[ i_p ] );
            println( "lp[ i_p + 1 ] : ", lp[ i_p + 1] );
            println( "p_loc         : ", p_loc );
            println( " i_p : ", i_p );
            println( " len : ", length( lp ) );
            println( prm[ i - 1 ] );
            println( prm[ i     ] );
            println( prm[ i + 1 ] );
#            println( prm );

        end
        #println( lp[ i_p ], "<=", p_loc );
        @assert( lp[ i_p ] <= p_loc  );
        while  ( i_q < length( lq ) )  &&  ( curr[ 2 ] >= lq[ i_q + 1 ] )
            i_q = i_q + 1;
        end

        #        println( "get_point P..." );
        pcurr,t_p = get_point( P, lp, i_p, p_loc );
        @assert( 0 <= t_p  &&  t_p <= 1.0 )
        #        println( "get_point Q ..." );
        get_point( Q, lq, i_q, curr[ 2 ] );
#        println( "get_point Q!!!!!!!!!! ..." );
        qcurr,t_q = get_point( Q, lq, i_q, curr[ 2 ] );

        if  ( ! ( 0.0 <= t_q <= 1.0 ) )
            println( "t_q: ", t_q );
            @assert( 0.0 <= t_q <= 1.0 );
        end
        #        println( "here? " );

        if  ( t_p == 0.0 ) ||  ( ( t_p == 1.0 )  &&  ( i_p == len_p ) )
            push!( pes, EventPoint( pcurr, i_p, PT_VERTEX, 0.0 ) );
        else
            push!( pes, EventPoint( pcurr, i_p, PT_ON_EDGE, t_p ) );
        end

        if  ( t_q == 0.0 ) ||  ( ( t_q == 1 )  &&  ( i_q == len_q ) )
            push!( qes, EventPoint( qcurr, i_q, PT_VERTEX, 0.0 ) );
        else
            push!( qes, EventPoint( qcurr, i_q, PT_ON_EDGE, t_q ) );
        end
        #       println( "here? B" );
    end

    push!( pes, EventPoint( last(P), cardin( P ), PT_VERTEX, 0.0 ) );
    push!( qes, EventPoint( last(Q), cardin( Q ), PT_VERTEX, 0.0 ) );

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

    return  pes, qes
end

########################################################################
"""
    Morphing_from_prm( prm, P, Q )

    Construct a morphing from a parameterization of the two polygons P
    and Q.

"""
function  Morphing_from_prm( prm::Polygon2F,
                             P::Polygon{N,T},
                             Q::Polygon{N,T} ) where {N,T}
    pes, qes = event_sequences_extract( prm, P, Q );
    return  Morphing_init( P, Q, pes, qes );
end


function   compute_offsets( P::Polygon{N,T},
                            pes::Vector{EventPoint{N,T}},
                            qes::Vector{EventPoint{N,T}},
                            )  where  {N,T}
    offs = zeros( Float64, cardin( P ) );

    for  i  in eachindex( pes )
        ev = pes[ i ];
        evx = qes[ i ];
        r = Dist( ev.p, evx.p );
        offs[ ev.i ] = max( offs[ ev.i ], r )
    end

    return  offs
end

"""
    Morphing_extract_offsets
"""
function  Morphing_extract_offsets( m::Morphing )
    p_offs = Vector{Float64}();
    q_offs = Vector{Float64}();


    p_offs = compute_offsets( m.P, m.pes, m.qes );
    q_offs = compute_offsets( m.Q, m.qes, m.pes );
    ; # P event sequence

    return  p_offs, q_offs
end





###########################################################################
# Compute the morphing src(v) -> trg(u): That is u(v(t))
###########################################################################

"""
    Morphing_combine

Gets two morphings u, v (i.e., two parameterizations) and combine them
into a single morphing u(v(.)).

For example, if u: γ → δ  and  v: δ → ξ, then the returned morphing is
u(v(⋅)): γ → ξ.

"""
function  Morphing_combine( u::Morphing{N,T}, v::Morphing{N,T} ) where {N,T}
    u_prm = Morphing_extract_prm( u );
    v_prm = Morphing_extract_prm( v );

    @assert( cardin( u_prm ) > 1 );
    @assert( cardin( v_prm ) > 1 );

    #println( "b5678" );
    prm = parameterization_combine( v_prm, u_prm );

    pes, qes = event_sequences_extract( prm, u.P, v.Q );

    m = Morphing_init( u.P, v.Q, pes, qes );

    return  m;
end


function  Morphing_adtw_price( m::Morphing{N,T} ) where  {N,T}
    P,Q = Morphing_as_polygons( m );

    len = cardin( P );
    price::Float64 = 0;
    for  i in 1: len-1
        price = price + segs_match_price( P[ i ], P[ i + 1 ],
                                          Q[ i ], Q[ i + 1 ] );
    end

    return  price;
end


Morphing2F = Morphing{2,Float64};

#
# End of file
##########################################################################
