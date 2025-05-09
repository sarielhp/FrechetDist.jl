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

include( "DistFunc.jl" );


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
struct EventPoint{N,T}
    p::Point{N,T}
    i::Int64;
    type::FPointType;
    t::Float64;
    penalty::Float64;
end

function  ev_update_p_t( ev::EventPoint{N,T}, p::Point{N,T}, t::Float64
                         ) where {N,T}
    return  EventPoint( p, ev.i, ev.type, t, ev.penalty );
end
function  ev_update_penalty( ev::EventPoint{N,T}, penalty 
                         ) where {N,T}
    return  EventPoint( ev.p, ev.i, ev.type, ev.t, penalty );
end
# leash = maximum length of an edge in the encoded matching.
#         (That is, the Frechet distance.)
"""
    Morphing

Encoding of a morphing (i.e., matching) between two polygonal cuves.

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

    f_is_monotone_init::Bool
    f_is_monotone::Bool

    lower_bound::T;
end

function   Morphing_init( P::Polygon{N,T}, Q::Polygon{N,T},
    pes::Vector{EventPoint{N,T}},
    qes::Vector{EventPoint{N,T}} ) where {N,T}

    @assert( length( pes ) == length( qes ) );
    m = Morphing( P, Q, pes, qes,
        0.0, 0, 0.0, 0.0, 0.0, false, false,
        0.0 # lower_bound
    );
    Morphing_recompute_leash( m );

    return  m;
end


function  Morphing_get_max_edges_err( m::Morphing{D,T} ) where {D,T}
    mx::Float64 = 0.0;

    io = jo = -1;

    P = m.P;
    Q = m.Q;

    l_p = cardin( P );
    l_q = cardin( Q );

    len = length( m.pes );
    for  ind in  1:len
        i = m.pes[ ind ].i;
        j = m.qes[ ind ].i;

        if  ( i== l_p )  ||  ( j == l_q )
            continue;
        end

        d = max( Dist( P[ i ], Q[ j ] ),
                 Dist( P[ i ], Q[ j + 1 ] ),
                 Dist( P[ i + 1 ], Q[ j ] ),
                 Dist( P[ i + 1 ], Q[ j + 1 ] ) );
        ld = iseg_iseg_dist( P[ i ], P[ i + 1 ], Q[ j ], Q[ j + 1 ] );
        delta = d - ld;
        #println( "delta : ", delta, " (", i, ", ", j, ")" );
        if  ( delta > mx )
            mx = delta;
            io, jo = i, j;
        end
    end

    return   io, jo, mx;
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
        # BUG: If should not have been there... We are being very
        # conservative about vertices that should not be simplified.

        # if ( ev.type == PT_VERTEX )
        vl[ ev.i ] = max( vl[ ev.i ], Dist( ev.p, evx.p ) );
        #end
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
    return  Morphing( P, Q, pes, qes, r, 0, 0.0, 0.0, 0.0, false, false, 0.0 );
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
            up = deepcopy( ep );
            push!( ns, up )
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
                p = convex_comb( P[ nep.i ], P[ nep.i + 1 ], t );
                nep = ev_update_p_t( nep, p, t );
            end
            push!( ns, nep );
            j = j + 1
        end

        nep = deepcopy( s[ j ] );
        if  ( t > nep.t )
            ell = Dist( P[ nep.i ], P[ nep.i + 1 ] );
            delta = max( delta, ( t - nep.t ) * ell );
            p = convex_comb( P[ nep.i ], P[ nep.i + 1 ],
                                 t );
            nep = ev_update_p_t( nep, p, t );
        end
        push!( ns, nep );

        i = j + 1;
    end

    @assert( length( ns ) == length( s ) );
    return  ns, delta;
end


##############################################################################
# Make the matching (which might be potentially not monotone) into a monotone
# one.
##############################################################################
function  events_seq_get_monotone_leash( P::Polygon{N,T},
    s::Vector{EventPoint{N,T}},
    s_alt::Vector{EventPoint{N,T}},
    leash::T ) where {N,T}

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
        t = s[ i ].t
        while  ( ( j < len )
                 &&  ( s[ j + 1 ].type == PT_ON_EDGE )
                 &&  ( s[ j + 1 ].i == loc ) )
            nep = s[ j ];
            t = max( t, nep.t );
            if  ( t > nep.t )
                ell = Dist( P[ nep.i ], P[ nep.i + 1 ] );
                new_p = convex_comb( P[ nep.i ], P[ nep.i + 1 ],
                                     t  );
                leash = max( leash, Dist( new_p, s_alt[ j ].p ) );
            end
            j = j + 1
        end

        nep = s[ j ];
        if  ( t > nep.t )
            ell = Dist( P[ nep.i ], P[ nep.i + 1 ] );
            new_p = convex_comb( P[ nep.i ], P[ nep.i + 1 ], t );
            leash = max( leash, Dist( new_p, s_alt[ j ].p ) );
        end
        i = j + 1;
    end

    return  leash;
end


function  Morphing_recompute_leash( m::Morphing{N,T} ) where  {N,T}
    @assert( length( m.pes ) == length( m.qes ) );
    r = 0;
    for  i  in eachindex( m.pes )
        ell = (Dist( m.pes[ i ].p, m.qes[ i ].p )
                        - m.pes[ i ].penalty - m.pes[ i ].penalty );
        if   ell > r
            r = ell;
        end
    end
    if  ( m.leash != r )   &&   ( m.leash != 0 )
        println( "ERROR old leash:" , m.leash );
        println( "ERROR new leash:" , r );
        exit( -1 );
    end
    m.leash = r;
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

    Morphing_verify_valid( m );

    #leash::Float64 = 0.0;
    for  i  in eachindex( m.pes )
        push!( P, m.pes[ i ].p );
        push!( Q, m.qes[ i ].p );
        d = Dist( m.pes[ i ].p, m.qes[ i ].p )
        #=if  ( d > leash )
            leash = max( leash, d );
        end=#
    end

    return  P, Q
end


"""
    Morphing_as_polygons_w_times

    Returns time for each pair of vertices in the morphing, between 0 and 1.
"""
function  Morphing_as_polygons_w_times( m::Morphing{N,T} ) where  {N,T}
    P, Q = Morphing_as_polygons( m );
    t = Vector{Float64}();

    len = polygon.total_length( P ) + polygon.total_length( Q );
    if  len == 0.0
        return  P, Q, zeros( length( P ) );
    end
    c_z::Int64 = 0;
    min_delta::Float64 = len;
    for  i  in  1:length(P) - 1
        delta = ( Dist( P[ i ], P[ i + 1] )
                  + Dist( Q[ i ], Q[ i + 1] ) );
        if  ( delta == 0.0 )
            c_z = c_z + 1;
            continue;
        end
        min_delta = min( min_delta, delta );
    end
    fake_delta = min_delta/8.0;
    len = len + fake_delta * c_z;

    push!( t, 0.0 );
    for  i  in  1:length(P) - 1
        delta = ( Dist( P[ i ], P[ i + 1] )
                  + Dist( Q[ i ], Q[ i + 1] ) );
        if  ( delta == 0.0 )
            delta = fake_delta;
        end
        next_t = min( last(t) + ( delta / len ), 1.0 );
        push!( t, next_t );
    end
    if last( t ) != 1.0
        pop!( t )
        push!( t, 1.0 );
    end
    return  P, Q, t
end


function  Morphing_as_function_w_times( m::Morphing{N,T} ) where  {N,T}
    P, Q = Morphing_as_polygons( m );
    t = Vector{Float64}();

    len = polygon.total_length( P );
    if  len == 0.0
        return  P, Q, zeros( length( P ) );
    end
    push!( t, 0.0 );
    c_z = 0;
    fake_delta = len / (4 + length( P ) );
    for  i  in  1:length(P) - 1
        delta = Dist( P[ i ], P[ i + 1] ) / len;
        if  ( delta == 0.0 )
            c_z = c_z + 1;
            len = len + fake_delta;
        end
    end
    for  i  in  1:length(P) - 1
        if  ( Dist( P[ i ], P[ i + 1] ) == 0.0 )
            delta = fake_delta / len;
        else
            delta = Dist( P[ i ], P[ i + 1] ) / len;
        end

        push!( t, min( last(t) + delta, 1.0 ) );
    end
    if last( t ) != 1.0
        pop!( t )
        push!( t, 1.0 );
    end
    return  P, Q, t
end


function  polygons_get_loc_at_time( P::Polygon{D,T},
                                    Q::Polygon{D,T},
                                    times::Vector{T},
                                    t::T
                                    ) where {D,T}
    ( t <= 0.0 ) &&  return first( P ), first( Q );
    ( t >= 1.0 ) &&  return last( P ), last( Q );

    pos = searchsortedfirst( times, t );
    if  (! ( 1 < pos <= length( P ) ))
        println( pos );
        println( times[ pos ] );
        println( t ) ;
        @assert( 1 < pos <= length( P )  );
    end

    #=
    println( "POS : ", pos );
    (pos > 1 )  &&  println( "t[pos-1]: ", times[ pos - 1 ] )
    println( "t[pos]: ", times[ pos ] )
    println( "qeury: ", t );
    =#
    prev = pos - 1;

    while  ( prev > 1 )  &&  (times[ prev - 1 ] == times[ prev ] )
        println( "GOING BACK!" );
        prev = prev - 1;
    end
    while  ( pos < length( times ) )  &&  (times[ pos ] == times[ pos + 1 ] )
        pos = pos + 1;
    end
    @assert( times[ prev ] <= t <= times[ pos ] );

    if  ( ( prev + 1 ) == pos )
        delta = ( t - times[ prev ] ) / (times[ pos ] - times[ prev ] );
        p = convex_comb( P[ prev ], P[ pos ], delta )
        q = convex_comb( Q[ prev ], Q[ pos ], delta )
        #println( "DDD: ", Dist( P[prev], Q[prev ] ), " PREV ", prev,
        #         "times: ", times[prev ] );
        return  p, q;
    end

    # The morphing had stopped at this point, and we should return the
    # maximum in this duration...
    @assert( false );


    return  p, q;
end

function   Morphing_sample_uniformly( m::Morphing{N,T}, n::Int64 ) where {N,T}
    P, Q, times = Morphing_as_polygons_w_times( m );

    m = zeros(n,2*N);
    for i in 0:n-1
        t = i/(n-1)
        p,q = polygons_get_loc_at_time( P, Q, times, t );
        m[i+1, :] = [ p... ; q... ];
#        m[ i, : ] = [ p... ; p... ];
    end

    return  m;
end

function   Morphing_is_monotone( m::Morphing{N,T} ) where {N,T}
    if  m.f_is_monotone_init
        return  m.f_is_monotone
    end
    m.f_is_monotone = ( events_seq_is_monotone( m.pes )
                      && events_seq_is_monotone( m.qes ) )
    m.f_is_monotone_init = true;

    return  m.f_is_monotone
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
        return  m; # deepcopy( m );
    end
    P = m.P;
    Q = m.Q;

    pes_new, da = events_seq_make_monotone( P, m.pes );
    qes_new, db = events_seq_make_monotone( Q, m.qes );


    m_out = Morphing_init( P, Q, pes_new, qes_new );

    m_out.monotone_err = max( da, db );

    return  m_out;
end


function  Morphing_monotone_leash( m::Morphing{N,T} ) where {N,T}
    if  Morphing_is_monotone( m )
        return  m.leash; # deepcopy( m );
    end

    leash = m.leash;
    leash = events_seq_get_monotone_leash( m.P, m.pes, m.qes, leash );
    leash = events_seq_get_monotone_leash( m.Q, m.qes, m.pes, leash );

    return  leash;
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

function  check_no_nan( V::Vector{EventPoint{N,T}} ) where {N,T}
    for  ev in V
        @assert( ! isnan( ev.p[1] ) )
        @assert( ! isnan( ev.p[2] ) )
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

    check_no_nan( m.pes );
    check_no_nan( m.qes );
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
    end at (m,n).  The polygonal curve either have positive slope edge,
    or vertical or horizontal edges. It can be thought of as a
    piecewise linear function from [0,m] to [0,n]. Here m and n are the
    lengths of the two given polygons of P and Q, respectively.

"""
function  Morphing_extract_prm( m::Morphing{N,T} )::cg.Polygon2F where {N,T}
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

    out = cg.Polygon2F();
    for  i in eachindex(pps)
        push!( out, npoint( pps[ i ], qps[ i ] ) )
    end
    return  out;
end


function  eval_pl_func_on_dim( p::Point{N,T}, q::Point{N,T}, val::Float64,
                               d::Int64 )  where {N,T}
    @assert( p[ d ] <= q[ d ] );
    t = (val - p[d]) / (q[d] - p[d]);

    return  convex_comb( p, q, t );
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



##############################################################3
# parameterization_combine
# Given two parameterization f, g, compute f(g( . ))
###############################################################
function   parameterization_combine( f::Polygon2F,
                                     g::Polygon2F )::Polygon2F
    if  ( ! fp_equal( last( g )[2], last( f )[ 1 ] ) )
        println( last( g )[2], " != ",  last( f )[ 1 ] );
        @assert( fp_equal( last( g )[2], last( f )[ 1 ] ) )
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
        f_equal_yf_yg::Bool = fp_equal( yf, yg )
        if  ( f_equal_yf_yg )  &&  ( idf == l_f )  &&  ( idg == l_g )
            push!( out, npoint( g[ idg ][1], f[ idf ][2] ) )
            break;
        end
        if  ( f_equal_yf_yg   &&  ( idf < l_f )
              &&  ( f[ idf + 1 ][1] == yf ) )
            push!( out, npoint( g[ idg ][1], f[ idf ][2] ) )
            idf = idf + 1
            continue;
        end
        if  ( f_equal_yf_yg   &&  ( idg < l_g )
              &&  ( fp_equal( g[ idg + 1 ][2], yg ) ) )
            push!( out, npoint( g[ idg ][1], f[ idf ][2] ) )
            idg = idg + 1
            continue;
        end
        if  f_equal_yf_yg
            push!( out, npoint( g[ idg ][1], f[ idf ][2] ) )
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
            push!( out, npoint( xg, f[ idf ][ 2 ] ) )
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
            push!( out, npoint( g[ idg ][ 1 ], zf ) )
            idg = min( idg + 1, l_g );
            continue;
        end
        @assert( false );
    end

    # We monotononize the output because of floating point error...
    # Should really just get rid of tiny flatting point jitters.
    curr::Point2F = deepcopy( out[ 1 ] );
    for  i in 2:cardin( out )
        curr = Point_max( out[ i ], curr );
        out[ i ] = curr;
        #=for  j in 1:2
            curr[ j ] = max( curr[ j ], out[ i ][ j ] )
            out[ i ][ j ] = curr[ j ];
        end=#
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
    p = convex_comb( P[ i ], P[i + 1], t );
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
        curr = prm[ i ];
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

        #t_p::Float64;
        pcurr,t_p = get_point( P, lp, i_p, p_loc );
        @assert( 0 <= t_p  &&  t_p <= 1.0 )
        #@assert( ! isNaN(pcurr) );
        #        println( "get_point Q ..." );
        #get_point( Q, lq, i_q, curr[ 2 ] );
        #        println( "get_point Q!!!!!!!!!! ..." );
        qcurr,t_q = get_point( Q, lq, i_q, curr[ 2 ] );
        #@assert( ! isNaN( qcurr ) );

        if  ( ! ( 0.0 <= t_q <= 1.0 ) )
            println( "t_q: ", t_q );
            @assert( 0.0 <= t_q <= 1.0 );
        end
        #        println( "here? " );

        if  ( t_p == 0.0 ) ||  ( ( t_p == 1.0 )  &&  ( i_p == len_p ) )
            push!( pes, EventPoint( deepcopy( pcurr), i_p, PT_VERTEX, 0.0, 0.0 ) );
        else
            push!( pes, EventPoint( deepcopy(pcurr), i_p, PT_ON_EDGE, T(t_p),
                                    0.0 ) );
        end

        #println( pcurr );
        #println( qcurr );
        if  ( t_q == 0.0 ) ||  ( ( t_q == 1 )  &&  ( i_q == len_q ) )
            push!( qes, EventPoint( deepcopy(qcurr), i_q, PT_VERTEX, 0.0, 0.0 ) );
        else
            push!( qes, EventPoint( deepcopy(qcurr), i_q, PT_ON_EDGE, T(t_q), 0.0 ) );
        end
        #       println( "here? B" );
    end

    @assert( !point.isNaN( last( P ) ) );
    @assert( !point.isNaN( last( Q ) ) );
    
    push!( pes, EventPoint( deepcopy(last(P)), cardin( P ), PT_VERTEX, 0.0,
                            0.0 ) );
    push!( qes, EventPoint( deepcopy(last(Q)), cardin( Q ), PT_VERTEX, 0.0,
                            0.0 ) );

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

function  Morphing_SweepDist_approx_price( m::Morphing{N,T} ) where  {N,T}
    P,Q = Morphing_as_polygons( m );

    @assert( length( P ) == length( Q ) );
    len = cardin( P );
    price::Float64 = 0;
    for  i in 1:(len-1)
        price = price + segs_match_price( P[ i ], P[ i + 1 ],
                                          Q[ i ], Q[ i + 1 ] );
    end

    return  price;
end


function  Morphing_SweepDist_price( m::Morphing{N,T} ) where  {N,T}
    Morphing_verify_valid( m );

    P,Q = Morphing_as_polygons( m );

    @assert( length( P ) == length( Q ) );
    len = cardin( P );
    price::Float64 = 0.0;
    for  i in 1:(len-1)
        delta::Float64 = SweepDist_segs( P[ i ], P[ i + 1 ],
            Q[ i ], Q[ i + 1 ] )
        price = price + delta;
    end

    return  price;
end

function  morphing_profile( m::Morphing{D,T}, n::Int64 ) where  {D,T}
    f_debug::Bool = false;

    if  ( f_debug )
        Morphing_recompute_leash( m );
        println( "LEASH recomputed: ", m.leash );
    end
    PB,QB,times = Morphing_as_polygons_w_times( m );
    #println( times );
    if  ( f_debug )
        for  i  in 1:length( times ) - 1
            if  ( ! ( times[ i ] <= times[ i + 1 ] ) )
                println( i );
                println( times[ i ], ", ", times[ i + 1 ] );
                @assert( ( times[ i ] <= times[ i + 1 ] ) );
            end
        end;
    end;
    #println( length( times ) );
    lP = total_length( PB );
    X = range(0, lP, length=n)

    Y = Vector{Float64}();
    for x ∈ X
        pos = x / lP;
        #println( " pos : ", pos );
        p, q = polygons_get_loc_at_time( PB, QB, times, pos );
        push!( Y, Dist(p,q ) );
    end
    f_debug  &&  println( maximum( Y ) );
    return  X, Y;
end


Morphing2F = Morphing{2,Float64};

#
# End of file
##########################################################################
