# Originally contributed by S. Har-Peled
# under MIT License
####
#---------------------------------------------------------
# Basic computational geometry library
#
# Points, segments, etc
#---------------------------------------------------------

module cg



#using PrecompileTools
using Parameters
using StaticArrays

#include
include( "point.jl" );
using .point

include( "polygon.jl" );
using .polygon

using LinearAlgebra
using DelimitedFiles
using .point

#using .Point;
#using .Polygon;
#using LoopVectorization
#using SIMD




######################################################################
######################################################################
# Bounding box
"""
    BBox

Axis parallel bounding box.
"""

@with_kw mutable struct BBox{D,T}
    f_init::Bool = false
    mini::MVector{D,T} = zeros( T, D );
    maxi::MVector{D,T} = zeros( T, D );
end

function  BBox_width( bb::BBox{D,T}, dim::Int64 = 1 ) where {D,T}
    return  bb.maxi[ dim ] - bb.mini[ dim ];
end
function  BBox_height( bb::BBox{D,T}, dim::Int64 = 2 ) where {D,T}
    return  bb.maxi[ dim ] - bb.mini[ dim ];
end

function  BBox_print( bb::BBox{D,T}) where {D,T}
    for  i in 1:D
        print( " [[" );
        print( bb.mini[ i ] );
        print( ".." )
        print(  bb.maxi[ i ] );
        print(  "]] " );
    end
end


function  BBox_get_dim( bb, d )
    return  ( bb.mini[ d ], bb.maxi[ d ] )
end

function  BBox_middle( bb::BBox{D,T}, d ) where  {D,T}
    return   (bb.mini[ d ] + bb.maxi[ d ]) / 2.0;
end
function  BBox_min( bb::BBox{D,T}, d ) where  {D,T}
    return  bb.mini[ d ];
end
function  BBox_max( bb::BBox{D,T}, d ) where  {D,T}
    return  bb.maxi[ d ];
end

function  BBox_init( bb::BBox{D,T}, p, q ) where  {D,T}
    f_init = true;
    for  i in 1:D
        bb.mini[ i ] = min( p[i ], q[ i ] );
        bb.maxi[ i ] = max( p[i ], q[ i ] );
    end
end


function  BBox_diam( bb::BBox{D,T} )  where  {D,T}
    ( ! bb.f_init )  &&  return  zero( T );

    sum = zero( T );
    for  i in 1:D
        sum += (bb.maxi[ i ] - bb.mini[ i ])^2;
    end

    return  sqrt( sum );
end

function  BBox_dist_in_dim( b::BBox{D,T}, c::BBox{D,T}, i::Int64
                             )  where  {D,T}

    ( b.mini[ i ] > c.maxi[ i ] )  &&  return  ( b.mini[ i ] - c.maxi[ i ] );
    ( b.maxi[ i ] < c.mini[ i ] )  &&  return  ( c.mini[ i ] - b.maxi[ i ] );

    return  zero( T );
end

"""
    BBox_extent_in_dim

    The maximum distance between two points in this dimension, that
belongs to the two bounding boxes.
"""
function  BBox_max_dist_in_dim( b::BBox{D,T}, c::BBox{D,T}, i::Int64
                               )  where  {D,T}

    ( b.mini[ i ] >= c.maxi[ i ] )  &&  return  ( b.maxi[ i ] - c.mini[ i ] );
    ( b.maxi[ i ] <= c.mini[ i ] )  &&  return  ( c.maxi[ i ] - b.mini[ i ] );

    # Boxes must intersect.... 
    return  max( abs( b.mini[ i ] - c.mini[ i ] ),
                 abs( b.maxi[ i ] - c.maxi[ i ] ),
                 abs( b.mini[ i ] - c.maxi[ i ] ),
                 abs( b.maxi[ i ] - c.mini[ i ] ) );
end

function  BBox_dist( b::BBox{D,T}, c::BBox{D,T} )  where  {D,T}
    ( ! f_init )  &&  return  zero( T );

    sum = zero( T );
    for  i in 1:D
        sum += ( BBox_dist_in_dim( b, c, i ) )^2;
    end

    return  sqrt( sum );
end


function  BBox_max_dist( b::BBox{D,T}, c::BBox{D,T} )  where  {D,T}
    ( ! f_init )  &&  return  zero( T );

    sum = zero( T );
    for  i in 1:D
        sum += ( BBox_max_dist_in_dim( b, c, i ) )^2;
    end

    return  sqrt( sum );
end


function  BBox_bottom_left( bb::BBox{D,T} ) where  {D,T}
    return  Point{D,T}( bb.mini );
end

function  BBox_top_right( bb::BBox{D,T} ) where  {D,T}
    return  Point{D,T}( bb.maxi );
end

function  BBox_expand( bb::BBox{D,T}, factor ) where  {D,T}
    mid::MVector{D,T} = (bb.maxi + bb.mini) / 2.0;
    diff::MVector{D,T} = (bb.maxi - bb.mini) * (factor/2.0);

#    margin::MVector{D,T} = (bb.maxi - bb.mini) * factor;
    bb.mini = mid - diff;
    bb.maxi = mid + diff;
end


###############################################33
###############################################33
### Line type
###############################################33


"""
    Line in D dimensions.

    `p` is a point on the line and `u` is the direction vector (not
    necessarily normalized). Parametrised as \$p + ut\$

"""
struct Line{D,T}
    p::MVector{D,T}
    u::MVector{D,T}
end

###############################################33
###############################################33
### Segment type
###############################################33
#Point{D,T} = point.Point{D,T}

"""
    Segment

Specifies a *directed* segment by two endpoints.
"""
struct Segment{D,T}
    p::Point{D,T}
    q::Point{D,T}
end

function  Segment{D,T}() where{D,T}
    return   Segment{D,T}(Point{D,T}(), Point{D,T}() );
end

function  Segment_get_on( seg::Segment{D,T}, t::Float64 ) where{D,T}
    return convex_comb( seg.p, seg.q, t );
end



function Base.show(io::IO, s::Segment{D,T}) where {D,T}
    print( io," [[" );
    print( io, s.p );
    print( io, ".." )
    print( io, s.q );
    print( io, "]] " );
end


function  Segment_length( seg::Segment{D,T} ) where {D,T}
    return point.Dist( seg.p, seg.q );
end

function  Dist( s::Segment{D,T}, qr::Point{D,T} ) where {D,T}
    nn = Segment_nn_point( s, qr );
    return  point.Dist( nn, qr );
end

function  Dist( a::Segment{D,T}, b::Segment{D,T} ) where {D,T}
    return  iseg_iseg_dist( a.p, a.q, b.p, b.q );
end

function  Segment_get_convex_coef( s::Segment{D,T}, qr::Point{D,T} ) where{D,T}
    d = point.Dist( s.p, qr );
    len = Segment_length( s );
    if  len == 0.0
        return  0.0;
    end
    return  d/len;
end

# Get nearest point on segment...
"""
    iseg_nn_point

Returns the closest point to the segment induced by the first two
points, to the query point. By avoiding creating the segment iself, it
is hopeflly more efficient.
"""
function  iseg_nn_point( s_p::Point{D,T}, s_q::Point{D,T},
                                qr::Point{D,T} ) where {D,T}
    # v(t) = p*(1-t) + q * t
    #      = p + t*(q-p)
    # v( [0,1] ) = segment.
    # u(t) = v(t) - qr = (p-qr) + t * (q-p)
    # D(t) = ||u(t)||^2 = ||p-qr||^2 + 2t <p-qr,q-p> + t^2 ||q-p||^2
    # Settig: a =  ||q-p||^2  and  b =  2 <p-qr,q-p>
    #        So D(t) = a*t^2 + b t + c
    # The minimum distance is achived at
    #    t^* = -b/(2a).
    a = DistSq( s_p, s_q );
    b = 2.0* dot( sub(s_p, qr), sub(s_q, s_p) );
    t::T = -b /(2.0 * a);

    if  ( t < 0 )
        t = 0;
    end
    if  ( t > 1 )
        t = 1;
    end
    pon = convex_comb( s_p, s_q, t );

    lon = DistSq(qr, pon);
    lp = DistSq(qr, s_p);
    lq = DistSq(qr, s_q);

    if  (( lon < lp )  &&  ( lon < lq ) )
        return  pon;
    end
    if  ( lp < lq )
        return  s_p;
    else
        return  s_q;
    end
end


function  iseg_iseg_dist( a_p::Point{D,T}, a_q::Point{D,T},
    b_p::Point{D,T}, b_q::Point{D,T} ) where {D,T}
    # v(s,t) = a_p*(1-s) + a_q * s - b_p*(1-t) - b_q * t
    #      = (a_p - b_p)  +  s * (a_q-a_p)  +  t * (b_p-b_q)
    #      = v_1  +  s * v_2  +  t * v_3
    #v_1::Point{D,T};
    #v_2::Point{D,T};
    #v_3::Point{D,T};

    v_1 = sub(a_p, b_p);
    v_2 = sub(a_q, a_p);
    v_3 = sub(b_p, b_q);

    #=

    D(s,t)
    = <v_1 + s*v_2 + t *v_3, v_1 + s*v_2 + t *v_3>
    = ||v_1||^2  +  2 * s * <v_1, v_2>  +  2 * t *<v_1,v_3>
       + s^2 * ||v_2||^2 + 2*s*t*<v_2,v_3>  + t^2 ||v_3||^2
    =
    Need to solve the linear system:

    0 = D'_s(s,t) = 2*<v_1,v_2> + 2*s*||V_2||^2  +  2*t<v_2,v_3>
    0 = D'_t(s,t) = 2*<v_1,v_3> + 2*s*<v_2,v_3> + 2*t * ||v_3||^2

    or equivalently:

    -<v_1,v_2> = s * ||V_2||^2  +  t * <v_2,v_3>
    -<v_1,v_3> = s * <v_2,v_3>  +  t * ||v_3||^2

    =#

    # Coefficients
    c = [ (- cg.dot( v_1, v_2 )),   (-cg.dot( v_1, v_3 ) ) ];

    m = [ ( cg.dot(v_2, v_2))  ( cg.dot(v_2, v_3) );
          ( cg.dot(v_2, v_3))  ( cg.dot(v_3, v_3) ) ];

    rk = rank( m );
    @assert( rk > 0 );
    s::Float64 = 0.0;
    y::Float64 = 0.0;
    #println( "rk = ", rk );
    if  ( rk == 1 )
        # The minimum distance is realized by one of the endpoints.
        d = min(
            dist_iseg_nn_point( a_p, a_q, b_p ),
            dist_iseg_nn_point( a_p, a_q, b_q ),
            dist_iseg_nn_point( b_p, b_q, a_p ),
            dist_iseg_nn_point( b_p, b_q, a_q ) );
        return  d;
    end
    #=
        println( "--- a, b --------------------" );
        println( a_p, a_q );
        println( b_p, b_q );

        println( "--- v_1,2,3 --------------------" );

        println( v_1 );
        println( v_2 );
        println( v_3 );

        println( "\n\n--- m --------------------" );
        println( "rank m: ", rank( m ) );
        println( m );
        println( "\n\n--- c --------------------" );
        println( c );

        println( "m.size: ", size( m ) );
        println( "c.size: ", size( c ) );
        =#

    # Solve the system...
    b = m \ c;
    #println( "b.size: ", size( b ) );
    s = b[ 1 ];
    t = b[ 2 ];

    # Snap solution if needed to the [0,1.0] interval....
    s = max( min( s, 1.0 ), 0.0 )
    t = max( min( t, 1.0 ), 0.0 )

    d::Float64 =  Dist( convex_comb( a_p, a_q, s ),
        convex_comb( b_p, b_q, t ) );

    d = min( d, dist_iseg_nn_point( a_p, a_q, b_p ),
                dist_iseg_nn_point( a_p, a_q, b_q ),
                dist_iseg_nn_point( b_p, b_q, a_p ),
                dist_iseg_nn_point( b_p, b_q, a_q ) );


    da = Dist( a_p, b_p );
    if  ( da < d )  &&  ( abs( da - d ) > (1e-8 * (d + da ) ) )
        println( "da < d? " );
        println( "da: ", da );
        println( "d ", d );
        println( "s: ", s );
        println( "t: ", t );
        println( "b: ", b );
        exit( -1 );

    end
    db = Dist( a_q, b_q );
    if  ( db < d )  &&  ( abs( db - d ) > (1e-8 * (d + db ) ) )
        println( "db < d? " );
        println( "db: ", da );
        println( "d ", d );
    end

    #println( "d = ", d );

    return  d;
end



#############################################
# Get nearest point on segment...
#
# s_p - s_q : the segment
# qr : the query point.
#############################################
"""
    dist_iseg_nn_point

Returns the *distance* to the closest point lying on the segment induced by
the first two points, to the query point. By avoiding creating the
segment iself, it is hopeflly more efficient.
"""
function dist_iseg_nn_point( s_p::Point{D,T}, s_q::Point{D,T}, qr::Point{D,T}
)  where {D,T}
    # v(t) = p*(1-t) + q * t
    #      = p + t*(q-p)
    # v( [0,1] ) = segment.
    # u(t) = v(t) - qr = (p-qr) + t * (q-p)
    #
    # D(t) = ||u(t)||^2
    #      = ||p-qr||^2 + 2t <p-qr,q-p> + t^2 ||q-p||^2
    #
    # Settig: a =  ||q-p||^2  and  b =  2 <p-qr,q-p>
    #        So D(t) = a*t^2 + b t + c
    # The minimum distance is achived at
    #    t^* = -b/(2a).
    #dff = Point{D,T}()

    #sub_dst( s_p, qr, dff );
    dff = sub( s_p, qr );

    a = DistSq( s_p, s_q );
    b = 2.0* dot( dff, sub( s_q, s_p ) );

    c = dot( dff, dff );

    t = -b /(2.0 * a);

    if  ( t < 0 )
        t = 0;
    end
    if  ( t > 1 )
        t = 1;
    end
    val::Float64 = a*t^2 + b*t + c;
    if  ( val <= 1e-20 )
        return 0;
    end
    return sqrt( val );
end


"""
    Segment_nn_point

Returns the closest point on the segment `s` to the query point `qr`.

"""
function  Segment_nn_point( s::Segment{D,T}, qr::Point{D,T} ) where {D,T}
    return   iseg_nn_point( s.p, s.q, qr );
end


#######################################################################
# Check if the plane bisector of p and q intersect seg, and if so where...
#f_on,t,p
#######################################################################
"""
    Segment_get_bisection_point -> Bool, Real, Point

    Computes the intersection point of the segment `seg` with the
bisector plane between `p` and `q`.

# Returns

The first argument returns whether the segment intersects the
bisector, the pramaterized location (tm), and the intersection piont
itself.

"""
function  Segment_get_bisection_point( seg::Segment{D,T}, p::Point{D,T},
    q::Point{D,T}  ) where {D,T}
    # Consider the segment going through p and q, and its middle point (mid).
    dir = sub( q, p );
    mid = ( q + p ) /2;
    pos = dot( dir, mid );

    # The line induced by seg, we are interested in the t, when
    # its dot proct with dir is equal to pos.

    # p::Point{D,T}
    # q::Point{D,T}
    vec = seg.q - seg.p;

    # seg.p + tm * vec:  The line under consideration
    # So:
    #     dot( dir, seg.p ) + t * dot( dir, vec ) = pos
    tm::Float64 = ( pos - dot( dir, seg.p ) ) / dot( dir, vec ) ;

    f_on::Bool = (0.0 <= tm <= 1.0);
    out = Segment_get_on( seg, tm );

#    println( "" );
#    println( "DDDD 1:", Dist( p, out ), "  D2: ", Dist( q, out ) );
#    println( "" );

    return  f_on, tm, out;
end




function Base.show(io::IO, p::Point{D,T}) where {D,T}
    print( io, "(" );
    for i in 1:D
        print( io, p[i] );
        if   i < D
            print( io, ", " );
        end
    end
    print( io, ")" );
end


function Base.show(io::IO, poly::Polygon{D,T}) where {D,T}
    f_iter::Bool = false;
    for p in Points( poly)
        if  f_iter
            print( io, "-" );
        else
            f_iter = true;
        end
        show(io, p );
    end
end

# get_max_pairwise_distance
function  DistInfty( P::Polygon{D,T},
                     Q::Polygon{D,T} )  where {D,T}

    n_p = cardin( P );
    n_q = cardin( Q );

    @assert  n_p == n_q

    dist = 0;
    for  i  in  1:n_p
        ell = Dist( P[i], Q[i] );
        if  ( ell > dist )
            dist = ell
        end
    end

    return  dist
end







#function  vert( p::Polygon{D,T}, i::Int64 ) where {D,T}
#    return  p.pnts[ i ];
#end



# Returns a polygon where the points are sampled uniformly along the polygon
#=
function  Polygon_sample_uniformly_old( poly::Polygon{D,T}, n::Int64 ) where {D,T}
    prefix_lens = Polygon_prefix_lengths( poly );
    len = last( prefix_lens );

    delta = len / (n-1);

    new_poly::Polygon{D,T} = Polygon{D,T}( Point{D,T}[] );
    push!( new_poly.pnts, first( poly.pnts ) );

    for  i  in  2:n-1
        t = delta * (i-1);
        p = Polygon_get_point_on( poly, prefix_lens, t );
        push!( new_poly.pnts, p );
    end
    push!( new_poly.pnts, last( poly.pnts ) );

    return  new_poly;
end
=#


function  Polygon_spine( P::Polygon{D,T} )  where {D,T}
    pout = Polygon{D,T}();
    if  ( cardin( P ) <= 0 )
        return  pout;
    end
    if  ( cardin( P ) <= 1 )
        Polygon_push( pout, P[ 1 ] );
        Polygon_push( pout, P[ 1 ] );
        return  pout;
    end
    Polygon_push( pout, P[ 1 ] );
    Polygon_push( pout, last( P ) );

    return  pout
end

"""
    Polygon_sample

Randomly sample each vertex of input polygon, with probability prob.
"""
function  Polygon_sample( P::Polygon{D,T}, prob::Float64 ) where  {D,T}
    out = Polygon{D,T}();
    for  p in P
        if  rand() <= prob
            push!( out, p )
        end
    end

    return  out;
end





distance(p1::Point{D,T}, p2::Point{D,T}) where {D,T} = norm( sub( p1, p2) )

function distance(y::Point{D,T}, l::Line{D,T}) where {D,T}
    p, u = l.p, l.u

    t = dot( (y - p), u ) / dot( u, u )
    x = Point(p + t*u)

    return Dist(x, y)
end


#######################################################################3

function  BBox_bound(  bb::BBox{D,T}, pnt::Point{D,T} )  where  {D,T}
    if  ! bb.f_init
#        println( "BB INIT!" );
        bb.f_init = true;
        bb.mini = deepcopy( pnt );
        bb.maxi = deepcopy( pnt );
        return
    end

    for  i in 1:D
        if  pnt[ i ] < bb.mini[ i ]
            bb.mini[ i ] = pnt[ i ];
        end
        if  pnt[ i ] > bb.maxi[ i ]
            bb.maxi[ i ] = pnt[ i ];
        end
    end
end

function  BBox_bound(  bb::BBox{D,T}, P::Polygon{D,T} )  where  {D,T}

    nv = cardin( P );
    for  i in 1:nv
        BBox_bound( bb, P[ i ] )
#        println( "i:", i, "   " );
#        BBox_print( bb );
#        println( "\n" );
    end
end

"""
    BBox_init

"""
function  BBox_init( P::Polygon{D,T},
                     range::UnitRange{Int64} ) where  {D,T}
    bb = BBox{D,T}()
    for  i in range
        BBox_bound( bb, P[ i ] );
    end

    return  bb;
end

function  BBox_bound(  bb::BBox{D,T},
                       list::Vector{Polygon{D,T}} )  where  {D,T}
    for  x in list
        BBox_bound( bb, x )
    end
end



#####################################################################
# Predefined useful types...



Segment2F = Segment{2,Float64};

VecPoint2I = Vector{Point{2,Int64}};
VecPolygon2F = Vector{Polygon{2,Float64}};
BBox2F = BBox{2,Float64};

#####################################################################



function  Polygon_random( D,T,  n::Int64 )
    #x::T = zero(T);
    #x = x + x;

    P = Polygon{D,T}();
    for  i in 1:n
        p::Point{D,T} = Point_random( D, T );
        push!( P, p  );
    end

    return P;
end

function  Polygon2F_random( n::Int64 )
    P = Polygon2F();

    for  i in 1:n
        v = rand( 2 );
        push!( P, point( v[1], v[2] ) );
    end
    return P;
end


"""
    Polygon_read_plt_orig_file

    Reads a .plt file into a polygon (2d floating point).
"""
function  Polygon_read_plt_orig_file( filename, dchar = "," )
    P::Polygon2F = Polygon2F();

    a = readdlm( filename );
    d_a::Int64 = size(a,1)
    #println( typeof( a ) );
    #println( size( a ) )
    # println( size(a,1) );
    # println( size(a,2) );

    for  r  in 7:size(a,1)
        line = a[r,1]
        #println( "line: [", line, "]" );

        parts::Vector{SubString{String}} = split.( line, dchar );
        x = parse( Float64, parts[ 1 ] );
        y = parse( Float64, parts[ 2 ] );
        push_smart!( P, point( x, y ) );
    end

    return  P
end

function get_numbers_in_line( s, ch )
    pieces = split(s, ch, keepempty=false)
    map(pieces) do piece
        #println( "Piece: [", piece, "]  ch[", ch, "]" );
        parse(Float64, piece)
    end
end


function  Polygon_read_file( filename, dchar = "," )
    P::Polygon2F = Polygon2F();

    lines = readlines( filename )
    for  i  in 1:length( lines )
        line = lines[ i ];

        # the sigspatial file have a header line... Just ignore it.
        if  ( line == "x y k tid" )
            #println( "CONTINUE!" );
            continue;
        end
        local pieces;
        if  occursin( " ", line )
            pieces = get_numbers_in_line( line, ' ' );
        else
            pieces = get_numbers_in_line( line, ',' );
        end
        @assert( length( pieces ) >= 2 );
        push_smart!( P, point( pieces[ 1 ], pieces[ 2 ] ) );
    end

    if  ( cardin( P ) <= 1 )
        Polygon_push( P, first( P ) );
    end
    return  P
end


function  Polygon_read_txt_file( filename, dchar = "," )
    P::Polygon2F = Polygon2F();

    println( "here!" );
    a = readdlm( filename );
    d_a::Int64 = size(a,1)
    println( typeof( a ) );
    println( size( a ) )
     println( size(a,1) );
     println( size(a,2) );
    exit(-1);

    return  P
end

Centroid( P ) = sum( P ) / length( P );


"""
    segs_match_price

The price of matching the edge p_a-p_b to the edge q_a-q_b.
"""
function  segs_match_price( p_a::Point{N,T},
                      p_b::Point{N,T},
                      q_a::Point{N,T},
                      q_b::Point{N,T} ) where {N,T}

    l_a = Dist( p_a, q_a );
    l_b = Dist( p_b, q_b );

    l_p = Dist( p_a, p_b );
    l_q = Dist( q_a, q_b );

    l_avg = ( l_a + l_b ) / 2.0;

    return   ( l_p + l_q ) * l_avg;
end

#####################################################################

export Segment
#export Polygon, Polygon2I, Point2I
export BBox
#export npoint
export BBox2F, Segment2F
#, Polygon2F, Point2F


########################################################
# Basic operations on points

export  add, sub, mult, norm

export   is_left_turn
export   is_right_turn


#
export  BBox_init, BBox_bound, BBox_expand, BBox_print, BBox_width
export  cardin, VecPolygon2F
export  BBox_bottom_left, BBox_top_right
export  BBox_min, BBox_max, BBox_middle;
export  BBox_diam;
export  BBox_dist, BBox_max_dist

export  iseg_nn_point

export  iseg_iseg_dist

export  Segment_nn_point
export  Segment_get_on, Segment_get_convex_coef
export  Segment_length
export  Segment_get_bisection_point

export  Points;
export  Centroid;

export  VecPnts_as_matrix

export  dist_iseg_nn_point

export  convex_comb

export  segs_match_price

end
