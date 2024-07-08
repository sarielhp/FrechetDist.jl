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
using LinearAlgebra
using DelimitedFiles
#using LoopVectorization

###################################################################
### Point type

"""
    Point

    Point in D dimensions. Implemented currently as a struct with
    StaticArray for values. It is templated, with `N` for dimension,
    and `T` for underlying type.

"""
# POINT
# mutable struct Point{D,T}
#    x::MVector{D,T}
# end

Point{D,T} = MVector{D,T};

function Point{D,T}() where  {D,T}
#    x::MVector{D,T} = zeros( T, D );
#    return  Point{D,T}( x );
    return  Point{D,T}( zeros( T, D ) );
end

#= POINT 
function  Base.getindex(p::Point{D,T}, i::Int) where {D,T}
   return   p.x[ i ];
end
=#

#= POINT
function  Base.setindex!(p::Point{D,T}, v, i::Int) where {D,T}
    p.x[ i ] = v;
end
=#

function  norm(p::MVector{D,T} ) where {D,T}
    sum = 0;
    for i in 1:D
        sum +=  p[i]^2;
    end
    return  sqrt( sum );
end

#=
function  norm(p::Point{D,T} ) where {D,T}
    return  norm( p.x );
end
=#

function  mult(z::T, p::Point{D,T}) where  {D,T}
    return  z * p;
    ### POINT    return  Point( z * p.x );
end
function  mult(p::Point{D,T}, z::T) where  {D,T}
    return  p * z;
    # POINT    return  Point( z * p.x );
end

#function  Base.:/( p::Point{D,T}, z ) where  {D,T}
function  pnt_div( p::Point{D,T}, z::T ) where  {D,T}
    return  p / z;
    # POINT    return  Point( (one(T) / z) * p.x );
end

function  normalize(p::Point{D,T} ) where {D,T}
    x = norm( p );
    if  x == 0
        return  p;
    end

    p = p*(1/x)
    return  p;
end



#function  Base.:-(p::Point{D,T}, q::Point{D,T}) where  {D,T}
function  sub(p::Point{D,T}, q::Point{D,T}) where  {D,T}
    #= POINT u = p.x - q.x;
    u = p - q;
    return  Point( u ); =#
    return  p - q;
end


#function  Base.:+(p::Point{D,T}, q::Point{D,T}) where  {D,T}
function  add(p::Point{D,T}, q::Point{D,T}) where  {D,T}
    #= POINT u = p.x + q.x;
    return  Point( u );
    =#
    return  p + q;
end

#function  Base.:*(p::Point{D,T}, c::T) where  {D,T}
#    u = p.x * c;
#    return  Point( u );
#end
#function  Base.:*(p::Point{D,T}, c::S) where  {D,T,S}
#    u = p.x * c;
#    return  Point( u );
#end

#function  dot( p::Point{D,T}, q::Point{D,T}) where  {D,T}
#    return  LinearAlgebra.dot( p, q );
#end
function  dot( p::Point{D,T}, q::Point{D,T}) where  {D,T}
    s = zero( T );
    for i in 1:D
        s += p[ i ] * q[ i ]; 
    end
    return  s;
    #return  LinearAlgebra.dot( p, q );
end

# Define standard lexicographical ordering on points
function Base.isless( p::Point{D,T}, q::Point{D,T} ) where {D,T}
    for i in 1:D
        if  p[ i ] < q[ i ]
            return  true;
        else
            if  p[ i ] > q[ i ]
                return  false;
            end
        end
    end
    return  false;
end

function  DistSq(p::Point{D,T}, q::Point{D,T}) where {D,T}
    sum = 0.0;
    for  i in 1:D
        sum += ( p[i] - q[i] )^2
    end

    return  sum  #norm( p1.x - p2.x );
end
function  Dist(p::Point{D,T}, q::Point{D,T}) where {D,T}
    sum = 0.0;
    for  i in 1:D
        sum +=  ( p[i] - q[i] )^2
    end

    return  sqrt( sum ) #norm( p1.x - p2.x );
end

function  convex_comb( p::Point{D,T}, q::Point{D,T}, t::Float64 ) where{D,T}
    if  ( 0 <= t <= 0.000001 )
        return  p;
    end
    if  ( 0.99999 < t <= 1.0 )
        return  q;
    end

    s = 1.0 - t;
    o = Point{D,T}(undef)

    @inbounds for  i in 1:D
        o[ i ] = p[ i ] * s + q[ i] * t;
    end

    return  0;
#    return  add( mult( p, 1.0-t), mult( q, t ) );
end


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

function  BBox_init( bb::BBox{D,T}, p, q ) where  {D,T}
    f_init = true;
    for  i in 1:D
        bb.mini[ i ] = min( p[i ], q[ i ] );
        bb.maxi[ i ] = max( p[i ], q[ i ] );
    end
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


"""
    point( args... )

A flexible constructor for a point specified by the arguments. Thus
point( 2.0, 3.0, 4.0 ) defined the 3d point (2.0, 3.0, 4.0). Or
similarly, point( 2.0, 1.0 ) would create a 2d point.

"""
function point( args...)
    D=length(args);
    T=typeof( first( args ) )
    x::MVector{D,T} = MVector{D,T}(undef);

    for i in eachindex( x )
        x[ i ] = args[ i ];
    end

    p::Point{D,T} = Point{D,T}( x );

    return  p;
end


function  Point_random( D, T )
    x::MVector{D,T} = MVector{D,T}(undef);

    for i in eachindex( x )
        x[ i ] = rand();
    end

    p::Point{D,T} = Point{D,T}( x );

    return  p;
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
    return Dist( seg.p, seg.q );
end

function  Dist( s::Segment{D,T}, qr::Point{D,T} ) where {D,T}
    nn = Segment_nn_point( s, qr );
    return  Dist( nn, qr );
end

function  Dist( a::Segment{D,T}, b::Segment{D,T} ) where {D,T}
    return  iseg_iseg_dist( a.p, a.q, b.p, b.q );
end

function  Segment_get_convex_coef( s::Segment{D,T}, qr::Point{D,T} ) where{D,T}
    d = Dist( s.p, qr );
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
    v_1::Point{D,T};
    v_2::Point{D,T};
    v_3::Point{D,T};
    
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


###############################################33
###############################################33
### Polygon type
###############################################33

struct Polygon{D,T}
    pnts::Vector{Point{D,T}};
end


function Polygon{D,T}()  where {D,T}
    vec = Vector{Point{D,T}}();
    return  Polygon{D,T}( vec );
end



function  Polygon_translate!( P::Polygon{D,T}, v::Point{D,T} ) where {D,T}
    for  i in 1:length(P.pnts)
        # POINT p.x = p.x - v.x;
        P.pnts[i] = sub( P.pnts[i], v );
    end

    return P;
end


function  Polygon_split_single_edge( P::Polygon{D,T}, i::Int64 ) where {D, T}
    if  i == cardin( P )
        i = i - 1;
    end
    if  i == 0
        i = i + 1;
    end
#    println( "iiii: ", i );
    Q = Polygon{D,T}();
    for t in 1:i
        Polygon_push_smart( Q, P[ t ] );
    end
    Polygon_push_smart( Q, (P[ i ] + P[ i+1 ] ) / 2.0 );
    for t in i+1:cardin(P)
        Polygon_push_smart( Q, P[ t ] );
    end
    return  Q;
end


function  Polygon_convex_comb( P::Polygon{D,T},
                               Q::Polygon{D,T},
                               t::Float64 ) where {D,T}
    n = cardin( P );
    @assert( n == cardin( Q ) );
    O = Polygon{D,T}();

    for i  in 1:n
        push!( O, convex_comb( P[ i ], Q[ i ], t ) );
    end

    return  O;
end


function  Polygon_as_matrix( P::Polygon{D,T} ) where {D,T}
    m = zeros( T, D, length( P.pnts ) );
    for  i in 1:length( P.pnts )
        m[:,i] = P.pnts[ i ];
    end
    return  m;
end

"""
    Polygon_write_to_file

    Writes a plt file of the polygon.
"""
function  Polygon_write_to_file( P::Polygon{D,T}, filename ) where {D,T}
    m = Polygon_as_matrix( P );
    open( filename, "w") do io
        writedlm( io, m', ',' )
    end
end


function  VecPnts_as_matrix( v::Vector{Point{D,T}} ) where {D,T}
    m = zeros( T, D, length( v ) );
    for  i in 1:length( v )
        m[:,i] = v[ i ];
    end
    return  m;
end


function Base.last( poly::Polygon{D,T} ) where {D,T}
#    println( "last??? " );
    return  last( poly.pnts );
end

function Base.first( poly::Polygon{D,T} ) where {D,T}
#    println( "last??? " );
    return  first( poly.pnts );
end


function  Base.getindex(a::Polygon{D,T}, i::Int) where {D,T}
    return   a.pnts[ i ];
end

function  Base.setindex!(a::Polygon{D,T}, v, i::Int) where {D,T}
    a.pnts[ i ] = v;
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
    for p in poly.pnts
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




function  Polygon_push_smart( pout::Polygon{D,T}, p::Point{D,T} ) where  {D,T}
    if  ( cardin( pout ) == 0 )
        push!( pout, deepcopy( p ) );
        return  true;
    end
    if  ( Dist( last( pout ), deepcopy( p ) ) > 0.0 )
        push!( pout, deepcopy( p ) );
        return  true;
    end
    return  false
end



"""
    Polygon_split_edges

Output a polygon, where each edge is split in the middle by
introducing a vertex into it.
"""
function  Polygon_split_edges( P::Polygon{D,T} ) where {D,T}
    l = cardin( P );
    Q = Polygon{D,T}();
    if  l == 0
        return  Q;
    end;
    for  i in 1:l-1
        p::Point{D,T} = convex_comb( P[ i ], P[ i + 1 ], 0.5 );
        Polygon_push_smart( Q, P[ i ] );
        Polygon_push_smart( Q, p );
    end
    Polygon_push_smart( Q, P[ l ] );

    return  Q;
end


function  Polygon_push( pout::Polygon{D,T}, p::Point{D,T} ) where  {D,T}
    push!( pout, deepcopy( p ) );
end



function Polygon_move_to_origin( P::Polygon{D,T} ) where {D,T}
    pout::Polygon{D,T} = Polygon{D,T}();
    p = P[ 1 ]
    for q in P.pnts
        Polygon_push_smart( pout, q - p );
    end
    return  pout;
end

function  cardin( p::Polygon{D,T} )::Int64 where {D,T}
    return  length( p.pnts );
end


function  Polygon_simplify_ext( P::Polygon{D,T}, r::T ) where {D,T}
#::Tuple{Polygon{D,T},Vector{Int64}}
    pout = Polygon{D,T}();
    pindices = Vector{Int64}();

    len = cardin( P );
    if  ( len == 0 )
        return pout;
    end
    if  ( Polygon_push_smart( pout, P[1] ) )
        push!( pindices, 1 );
    end

    curr = P[1];
    for i in 2:len
        if  ( Dist( P[i], curr ) > r )
            curr = P[ i ];
            if  ( Polygon_push_smart( pout, P[i] ) )
                push!( pindices, i );
            end
        end
    end

    # We push the last vertex in event if it is a repetition, for
    # example if the curve is close, for example...
    push!( pout, P[ len ] )
    push!( pindices, len );

    return  pout, pindices;
end


function  Polygon_simplify( P::Polygon{D,T}, r ) where {D,T}
    return   Polygon_simplify_ext( P, r )[ 1 ];
end

####################################################################
# We specify for each vertex how much error it is willing to accept
# What we do here is pretty Conservative and silly. Should be enough
# as a first try.
####################################################################
function  Polygon_simplify_radii( P::Polygon{D,T}, r::Vector{T} ) where {D,T}
    @assert( cardin( P ) == length( r ) );

    pindices = Vector{Int64}();
    pout = Polygon{D,T}();

    len = cardin( P );
    if  ( len == 0 )
        return pout, pindices;
    end
    if  Polygon_push_smart( pout, P[1] )
        push!( pindices, 1 );
    end

    curr = P[1];
    curr_r = r[ 1 ];
    for i in 2:len
        curr_r = min( curr_r, r[ i ] );
        if  ( Dist( P[i], curr ) > curr_r )
            curr = P[ i ];
            if  ( i < len )
                curr_r = r[ i + 1 ];
            end
            if  Polygon_push_smart( pout, P[i] )
                push!( pindices, i );
            end
        end
    end
    if  Polygon_push_smart( pout, P[ len ] )
        push!( pindices, len );
    end
    #println( "STRLEN: ", length( pindices ) )
    if  ( length( pindices ) <= 1 )
        Polygon_push( pout, P[ len ] )
        push!( pindices, len );
    end
    
    return  pout, pindices;
end


#function  vert( p::Polygon{D,T}, i::Int64 ) where {D,T}
#    return  p.pnts[ i ];
#end

function Base.push!( c::Polygon{D,T}, p::Point{D,T}) where {D,T}
    push!(c.pnts, p);
end

function Base.pop!( c::Polygon{D,T}) where {D,T}
    pop!(c.pnts);
end

function  Polygon_length( poly::Polygon{D,T} ) where {D,T}
    len::Float64 = 0.0;
    for i in firstindex(poly.pnts) : lastindex(poly.pnts)-1
        len += Dist( poly[ i ], poly[ i + 1 ] );
    end

    return len;
end


function  Polygon_prefix_lengths( poly::Polygon{D,T}
)::Vector{Float64} where {D,T}
    v = Vector{Float64}();
    push!( v,  zero(Float64) );
    len::Float64 = 0.0;
    for i in firstindex(poly.pnts) : lastindex(poly.pnts)-1
        p = poly.pnts[ i ];
        q = poly.pnts[ i + 1 ];
        len += Dist( p, q );
        push!( v, Float64( len ) );
    end
    return v;
end


function  Polygon_edge_length( P::Polygon{D,T}, i::Int64 )  where {D,T}
    if  i < 0  ||  ( i >= cardin( P ) )
        return  0.0;
    end
    return  Dist( P[ i ], P[ i + 1 ] );
end


# poly: Poylgon
# prefix_len: Prefix lengths of the edges,
# t: Location of points to be computed, where t is in the range
#    [0, euclidean_length( poly ) ]
function  Polygon_get_point_on( poly, prefix_len, t )
    if  t <= 0
        return first( poly.pnts );
    end
    #println( prefix_len, " : ", t );
    if  t >= last( prefix_len )
        return  last( poly.pnts );
    end
    i = searchsortedfirst( prefix_len, t )
    #    println( "i= ", i );
    if  prefix_len[ i ] == t
        return  poly.pnts[ i ];
    end

#    if  ( i > 1 )
    @assert( prefix_len[ i - 1 ] <= t <= prefix_len[ i ] );
#    println( prefix_len[ i - 1 ], " <= ", t, " <= ", prefix_len[ i ] );
#    end
#    if  ( prefix_len
    #println( "prefix_len[ i ]: ",  prefix_len[ i ] );
#    println( "t              : ", t );
    # prefix_len[ i -1 ] < t <  prefix_len[ i ]
    delta = prefix_len[ i ] - prefix_len[ i - 1 ];  # length of segment
    x = (t - prefix_len[ i - 1 ]) / delta;

    return  convex_comb( poly.pnts[ i - 1 ], poly.pnts[ i ], x );
end


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


############################################################################
# Returns a polygon where the points are sampled uniformly along the polygon
# n is roughly the number of vertices one getds
############################################################################

function  Polygon_sample_uniformly( P::Polygon{D,T}, n::Int64 ) where {D,T}
    prefix_lens = Polygon_prefix_lengths( P );
    len = last( prefix_lens );

    delta = len / (n-1);

    new_P::Polygon{D,T} = Polygon{D,T}( );
    Polygon_push_smart( new_P, first( P.pnts ) );

    sz = cardin( P );

    for  i  in  1:(sz-1)
        s = Segment{D,T}( P[ i ], P[ i + 1 ] );

        ell = Segment_length( s );

        ns::Int64 = floor( Int64, ell / delta );
        for  j  in 1:ns
            Polygon_push_smart( new_P,
                                Segment_get_on( s, j / ( ns + 1 ) ) );
        end
        Polygon_push_smart( new_P, P[ i + 1 ] );
    end

    return  new_P;
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

function  BBox_bound(  bb::BBox{D,T},
                       list::Vector{Polygon{D,T}} )  where  {D,T}
    for  x in list
        BBox_bound( bb, x )
    end
end



#####################################################################
# Predefined useful types...

Point2I = Point{2,Int64};
Point2F = Point{2,Float64};

Polygon2I = Polygon{2,Int64};
Polygon2F = Polygon{2,Float64};

Segment2F = Segment{2,Float64};

VecPoint2I = Vector{Point{2,Int64}};
VecPolygon2F = Vector{Polygon{2,Float64}};
BBox2F = BBox{2,Float64};

#####################################################################



function  Polygon_random( D, T, n::Int64 )
    P = Polygon{D,T}();

    for  i in 1:n
        push!( P, Point_random( D, T ) );
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
        Polygon_push_smart( P, point( x, y ) );
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
        Polygon_push_smart( P, point( pieces[ 1 ], pieces[ 2 ] ) );
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
export Polygon, Polygon2I, Point, Point2I
export BBox
export point
export BBox2F, Segment2F, Polygon2F, Point2F

# Basic operations on points

export  add, sub, mult, norm

# 
export  BBox_init, BBox_bound, BBox_expand, BBox_print, BBox_width
export  cardin, Dist, DistSq, VecPolygon2F
export  BBox_bottom_left, BBox_top_right

export  iseg_nn_point

export  iseg_iseg_dist

export  Polygon_prefix_lengths

export  Segment_nn_point
export  Segment_get_on, Segment_get_convex_coef
export  Segment_length
export  Segment_get_bisection_point

export  Polygon_length, Polygon_move_to_origin
export  Polygon_sample_uniformly, Polygon_push_smart, Polygon_spine

export  Polygon_read_file
export  Polygon_read_plt_orig_file

export  Polygon_read_txt_file


export  Polygon_simplify, Polygon_push, DistInfty
export  Polygon_simplify_radii
export  Polygon_simplify_ext
export  Polygon_translate!
export  Polygon_get_point_on
export  Polygon_as_matrix
export  Polygon_write_to_file
export  Polygon_random
export  Polygon_convex_comb
export  Polygon_split_edges
export  Polygon_split_single_edge
export  Polygon_edge_length

export  VecPnts_as_matrix

export  dist_iseg_nn_point

export  convex_comb

export  segs_match_price

end
