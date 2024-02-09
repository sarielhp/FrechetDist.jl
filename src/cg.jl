module cg
############################################################
# Basic computational geometry entities
#
# Points, segments, etc
############################################################

#using PrecompileTools
using Parameters
using StaticArrays
using LinearAlgebra
using DelimitedFiles

###################################################################
### Point type

"""
    Point

    Point in N dimensions. Implemented currently as a struct with
    StaticArray for values. It is templated, with `N` for dimension,
    and `T` for underlying type.

"""
mutable struct Point{N,T}
    x::MVector{N,T}

end

function Point{N,T}() where  {N,T}
#    x::MVector{N,T} = zeros( T, N );
#    return  Point{N,T}( x );
    return  Point{N,T}( zeros( T, N ) );
end

function  Base.getindex(p::Point{N,T}, i::Int) where {N,T}
    return   p.x[ i ];
end

function  Base.setindex!(p::Point{N,T}, v, i::Int) where {N,T}
    p.x[ i ] = v;
end

function  norm(p::MVector{N,T} ) where {N,T}
    sum = 0;
    for i in 1:N
        sum = sum + p[i]^2;
    end
    return  sqrt( sum );
end

function  norm(p::Point{N,T} ) where {N,T}
    return  norm( p.x );
end

function  normalize(p::Point{N,T} ) where {N,T}
    x = norm( p );
    if  x == 0
        return  p;
    end

    p = p*(1/x)
    return  p;
end

function  Base.:*(z, p::Point{N,T}) where  {N,T}
    return  Point( z * p.x );
end
function  Base.:/( p::Point{N,T}, z ) where  {N,T}
    return  Point( (1.0/z) * p.x );
end


function  Base.:-(p::Point{N,T}, q::Point{N,T}) where  {N,T}
    u = p.x - q.x;
    return  Point( u );
end
function  Base.:+(p::Point{N,T}, q::Point{N,T}) where  {N,T}
    u = p.x + q.x;
    return  Point( u );
end

function  Base.:*(p::Point{N,T}, c::T) where  {N,T}
    u = p.x * c;
    return  Point( u );
end
function  Base.:*(p::Point{N,T}, c::S) where  {N,T,S}
    u = p.x * c;
    return  Point( u );
end

function  dot( p::Point{N,T}, q::Point{N,T}) where  {N,T}
    return  LinearAlgebra.dot( p.x, q.x );
end

# Define standard lexicographical ordering on points
function Base.isless( p::Point{N,T}, q::Point{N,T} ) where {N,T}
    for i in 1:N
        if  p.x[ i ] < q.x[ i ]
            return  true;
        else
            if  p.x[ i ] > q.x[ i ]
                return  false;
            end
        end
    end
    return  false;
end

function  DistSq(p::Point{N}, q::Point{N}) where {N}
    sum = 0.0;
    for  i in 1:N
        sum = sum + ( p.x[i] - q.x[i] )^2
    end

    return  sum  #norm( p1.x - p2.x );
end
function  Dist(p::Point{N}, q::Point{N}) where {N}
    sum = 0.0;
    for  i in 1:N
        sum = sum + ( p.x[i] - q.x[i] )^2
    end

    return  sqrt( sum ) #norm( p1.x - p2.x );
end



######################################################################
######################################################################
# Bounding box
@with_kw mutable struct BBox{N,T}
    f_init::Bool = false
    mini::MVector{N,T} = zeros( T, N );
    maxi::MVector{N,T} = zeros( T, N );
end

function  BBox_width( bb::BBox{N,T}, dim::Int64 = 1 ) where {N,T}
    return  bb.maxi[ dim ] - bb.mini[ dim ];
end
function  BBox_height( bb::BBox{N,T}, dim::Int64 = 2 ) where {N,T}
    return  bb.maxi[ dim ] - bb.mini[ dim ];
end

function  BBox_print( bb::BBox{N,T}) where {N,T}
    for  i in 1:N
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

function  BBox_init( bb::BBox{N,T}, p, q ) where  {N,T}
    f_init = true;
    for  i in 1:N
        bb.mini[ i ] = min( p.x[i ], q.x[ i ] );
        bb.maxi[ i ] = max( p.x[i ], q.x[ i ] );
    end
end

function  BBox_bottom_left( bb::BBox{N,T} ) where  {N,T}
    return  Point{N,T}( bb.mini );
end

function  BBox_top_right( bb::BBox{N,T} ) where  {N,T}
    return  Point{N,T}( bb.maxi );
end

function  BBox_expand( bb::BBox{N,T}, factor ) where  {N,T}
    mid::MVector{N,T} = (bb.maxi + bb.mini) / 2.0;
    diff::MVector{N,T} = (bb.maxi - bb.mini) * (factor/2.0);

#    margin::MVector{N,T} = (bb.maxi - bb.mini) * factor;
    bb.mini = mid - diff;
    bb.maxi = mid + diff;
end



function point( args...)
    N=length(args);
    T=typeof( first( args ) )
    x::MVector{N,T} = MVector{N,T}(undef);

    for i in eachindex( x )
        x[ i ] = args[ i ];
    end

    p::Point{N,T} = Point{N,T}( x );

    return  p;
end


###############################################33
###############################################33
### Line type
###############################################33


"""
    Line in N dimensions.

    `p` is a point on the line and `u` is the direction vector (not
    necessarily normalized). Parametrised as \$p + ut\$

"""
struct Line{N,T}
    p::MVector{N,T}
    u::MVector{N,T}
end

###############################################33
###############################################33
### Segment type
###############################################33

struct Segment{N,T}
    p::Point{N,T}
    q::Point{N,T}
end

function  Segment{N,T}() where{N,T}
    return   Segment{N,T}(Point{N,T}(), Point{N,T}() );
end

function  Segment_get_on( seg::Segment{N,T}, t::Float64 ) where{N,T}
    v::Point{N,T} = seg.p * (1.0-t) + seg.q * t;
end

function  convex_comb( p::Point{N,T}, q::Point{N,T}, t::Float64 ) where{N,T}
    if  ( 0 <= t <= 0.000001 )
        return  p;
    end
    if  ( 0.99999 < t <= 1.0 )
        return  q;
    end

    return  p * (1.0-t) + q * t;
end


function Base.show(io::IO, s::Segment{N,T}) where {N,T}
    print( io," [[" );
    print( io, s.p );
    print( io, ".." )
    print( io, s.q );
    print( io, "]] " );
end

function  Segment_length( seg::Segment{N,T} ) where {N,T}
    return Dist( seg.p, seg.q );
end

function  Dist( s::Segment{N,T}, qr::Point{N,T} ) where {N, T}
    nn = Segment_nn_point( s, qr );
    return  Dist( nn, qr );
end

function  Segment_get_convex_coef( s::Segment{N,T}, qr::Point{N,T} ) where{N,T}
    d = Dist( s.p, qr );
    len = Segment_length( s );
    return  d/len;
end

# Get nearest point on segment...
function  induced_seg_nn_point( s_p::Point{N,T}, s_q::Point{N,T},
                                qr::Point{N,T} ) where {N, T}
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
    b = 2.0* dot( s_p-qr, s_q - s_p );
    t = -b /(2.0 * a);

    if  ( t < 0 )
        t = 0;
    end
    if  ( t > 1 )
        t = 1;
    end
    pon = s_p*(1-t) + s.q * t;

    lon = DistSq(qr, pon);
    lp = DistSq(qr, s_p);
    lq = DistSq(qr, s.q);

    if  (( lon < lp )  &&  ( lon < lq ) )
        return  pon;
    end
    if  ( lp < lq )
        return  s_p;
    else
        return  s_q;
    end
end


#############################################
# Get nearest point on segment...
#
# s_p - s_q : the segment
# qr : the query point.
#############################################
function  dist_seg_nn_point( s_p::Point{N,T}, s_q::Point{N,T},
                             qr::Point{N,T} ) where {N, T}
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
    dff = s_p - qr;

    a = DistSq( s_p, s_q );
    b = 2.0* dot( dff, s_q - s_p );

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




# Get nearest point on segment...
function  Segment_nn_point( s::Segment{N,T}, qr::Point{N,T} ) where {N, T}
    # v(t) = p*(1-t) + q * t
    #      = p + t*(q-p)
    # v( [0,1] ) = segment.
    # u(t) = v(t) - qr = (p-qr) + t * (q-p)
    # D(t) = ||u(t)||^2 = ||p-qr||^2 + 2t <p-qr,q-p> + t^2 ||q-p||^2
    # Settig: a =  ||q-p||^2  and  b =  2 <p-qr,q-p>
    #        So D(t) = a*t^2 + b t + c
    # The minimum distance is achived at
    #    t^* = -b/(2a).
    a = DistSq( s.p, s.q );
    b = 2.0* dot( s.p-qr, s.q - s.p );
    t = -b /(2.0 * a);

    if  ( t < 0 )
        t = 0;
    end
    if  ( t > 1 )
        t = 1;
    end
    pon = s.p*(1-t) + s.q * t;

    lon = DistSq(qr, pon);
    lp = DistSq(qr, s.p);
    lq = DistSq(qr, s.q);

    if  (( lon < lp )  &&  ( lon < lq ) )
        return  pon;
    end
    if  ( lp <= lq )
        return  s.p;
    else
        return  s.q;
    end
end

#######################################################################
# Check if the plane bisector of p and q intersect seg, and if so where...
#f_on,t,p
#######################################################################
function  Segment_get_bisection_point( seg::Segment{N,T}, p, q  ) where {N,T}
    # Consider the segment going through p and q, and its middle point (mid).
    dir = q - p;
    mid = ( q + p ) /2;
    pos = dot( dir, mid );

    # The line induced by seg, we are interested in the t, when
    # its dot proct with dir is equal to pos.

    # p::Point{N,T}
    # q::Point{N,T}
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

struct Polygon{N,T}
    pnts::Vector{Point{N,T}};
end


function Polygon{N,T}()  where {N,T}
    vec = Vector{Point{N,T}}();
    return  Polygon{N,T}( vec );
end

function  Polygon_translate!( P::Polygon{N,T}, v::Point{N,T} ) where {N,T}
    for  p in P.pnts
        p.x = p.x - v.x;
    end

    return P;
end


function  Polygon_as_matrix( P::Polygon{N,T} ) where {N,T}
    m = zeros( T, N, length( P.pnts ) );
    for  i in 1:length( P.pnts )
        m[:,i] = P.pnts[ i ].x;
    end
    return  m;
end

function  Polygon_write_to_file( P::Polygon{N,T}, filename ) where {N,T}
    m = Polygon_as_matrix( P );
    open( filename, "w") do io
        writedlm( io, m )
    end
end


function  VecPnts_as_matrix( v::Vector{Point{N,T}} ) where {N,T}
    m = zeros( T, N, length( v ) );
    for  i in 1:length( v )
        m[:,i] = v[ i ].x;
    end
    return  m;
end


function Base.last( poly::Polygon{N,T} ) where {N,T}
#    println( "last??? " );
    return  last( poly.pnts );
end

function Base.first( poly::Polygon{N,T} ) where {N,T}
#    println( "last??? " );
    return  first( poly.pnts );
end


function  Base.getindex(a::Polygon{N,T}, i::Int) where {N,T}
    return   a.pnts[ i ];
end

function  Base.setindex!(a::Polygon{N,T}, v, i::Int) where {N,T}
    a.pnts[ i ] = v;
end


function Base.show(io::IO, p::Point{N,T}) where {N,T}
    print( io, "(" );
    for i in 1:N
        print( io, p.x[i] );
        if   i < N
            print( io, ", " );
        end
    end
    print( io, ")" );
end


function Base.show(io::IO, poly::Polygon{N,T}) where {N,T}
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
function  DistInfty( P::Polygon{N,T},
                     Q::Polygon{N,T} )  where {N,T}

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




function  Polygon_push_smart( pout::Polygon{N,T}, p::Point{N,T} ) where  {N,T}
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

function  Polygon_push( pout::Polygon{N,T}, p::Point{N,T} ) where  {N,T}
    push!( pout, deepcopy( p ) );
end



function Polygon_move_to_origin( P::Polygon{N,T} ) where {N,T}
    pout::Polygon{N,T} = Polygon{N,T}();
    p = P[ 1 ]
    for q in P.pnts
        Polygon_push_smart( pout, q - p );
    end
    return  pout;
end

function  cardin( p::Polygon{N,T}  ) where {N,T}
    return  length( p.pnts );
end


function  Polygon_simplify_ext( P::Polygon{N,T}, r ) where {N,T}
    pout = Polygon{N,T}();
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
    if ( Polygon_push_smart( pout, P[ len ] ) )
        push!( pindices, len );
    end

    return  pout, pindices;
end


function  Polygon_simplify( P::Polygon{N,T}, r ) where {N,T}
    return   Polygon_simplify_ext( P, r )[ 1 ];
end

####################################################################
# We specify for each vertex how much error it is willing to accept
# What we do here is pretty Conservative and silly. Should be enough
# as a first try.
####################################################################
function  Polygon_simplify_radii( P::Polygon{N,T}, r::Vector{T} ) where {N,T}
    @assert( cardin( P ) == length( r ) );

    pindices = Vector{Int64}();
    pout = Polygon{N,T}();

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

    return  pout, pindices;
end


#function  vert( p::Polygon{N,T}, i::Int64 ) where {N,T}
#    return  p.pnts[ i ];
#end

function Base.push!( c::Polygon{N,T}, p::Point{N,T}) where {N,T}
    push!(c.pnts, p);
end

function Base.pop!( c::Polygon{N,T}) where {N,T}
    pop!(c.pnts);
end

function  Polygon_length( poly::Polygon{N,T} ) where {N,T}
    len::Float64 = 0.0;
    for i in firstindex(poly.pnts) : lastindex(poly.pnts)-1
        len += Dist( poly[ i ], poly[ i + 1 ] );
    end

    return len;
end


function  Polygon_prefix_lengths( poly::Polygon{N,T} ) where {N,T}
    v::Vector{Float64} = [];
    push!( v, 0 );
    len = 0;
    for i in firstindex(poly.pnts) : lastindex(poly.pnts)-1
        p = poly.pnts[ i ];
        q = poly.pnts[ i + 1 ];
        len += Dist( p, q );
        push!( v, len );
    end
    return v;
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

    return  poly.pnts[ i - 1 ] * ( one(x) - x) +  poly.pnts[ i ] * x;
end


# Returns a polygon where the points are sampled uniformly along the polygon
#=
function  Polygon_sample_uniformly_old( poly::Polygon{N,T}, n::Int64 ) where {N,T}
    prefix_lens = Polygon_prefix_lengths( poly );
    len = last( prefix_lens );

    delta = len / (n-1);

    new_poly::Polygon{N,T} = Polygon{N,T}( Point{N,T}[] );
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


function  Polygon_spine( P::Polygon{N,T} )  where {N,T}
    pout = Polygon{N,T}();
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

function  Polygon_sample_uniformly( P::Polygon{N,T}, n::Int64 ) where {N,T}
    prefix_lens = Polygon_prefix_lengths( P );
    len = last( prefix_lens );

    delta = len / (n-1);

    new_P::Polygon{N,T} = Polygon{N,T}( Point{N,T}[] );
    Polygon_push_smart( new_P, first( P.pnts ) );

    sz = cardin( P );

    for  i  in  1:(sz-1)
        s = Segment{N,T}( P[ i ], P[ i + 1 ] );

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



distance(p1::Point{N}, p2::Point{N}) where {N} = norm(p1.x - p2.x)

function distance(y::Point{N}, l::Line{N}) where {N}
    p, u = l.p, l.u

    t = (y.x - p) ⋅ u / (u ⋅ u)
    x = Point(p + t*u)

    return Dist(x, y)
end


#######################################################################3

function  BBox_bound(  bb::BBox{N,T}, pnt::Point{N,T} )  where  {N,T}
    if  ! bb.f_init
#        println( "BB INIT!" );
        bb.f_init = true;
        bb.mini = deepcopy( pnt.x );
        bb.maxi = deepcopy( pnt.x );
        return
    end

    for  i in 1:N
        if  pnt.x[ i ] < bb.mini[ i ]
            bb.mini[ i ] = deepcopy(pnt.x[ i ]);
        end
        if  pnt.x[ i ] > bb.maxi[ i ]
            bb.maxi[ i ] = deepcopy( pnt.x[ i ] );
        end
    end
end

function  BBox_bound(  bb::BBox{N,T}, P::Polygon{N,T} )  where  {N,T}

    nv = cardin( P );
    for  i in 1:nv
        BBox_bound( bb, P[ i ] )
#        println( "i:", i, "   " );
#        BBox_print( bb );
#        println( "\n" );
    end
end

function  BBox_bound(  bb::BBox{N,T},
                       list::Vector{Polygon{N,T}} )  where  {N,T}
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



function  Polygon2F_random( n::Int64 )
    P = Polygon2F();

    for  i in 1:n
        v = rand( 2 );
        push!( P, point( v[1], v[2] ) );
    end
    return P;
end


"""
    Polygon_read_plt_file

    Reads a .plt file into a polygon (2d floating point).
"""
function  Polygon_read_plt_file( filename )
    P::Polygon2F = Polygon2F();

    a = readdlm( filename );
    for  r  in 7:size(a,1)
        line = a[r,1]
        parts::Vector{SubString{String}} = split.( line, "," );
        x = parse( Float64, parts[ 1 ] );
        y = parse( Float64, parts[ 2 ] );
        Polygon_push_smart( P, point( x, y ) );
    end

    return  P
end


#####################################################################

export Segment
export Polygon, Polygon2I, Point, Point2I
export BBox
export point
export BBox2F, Segment2F, Polygon2F, Point2F
#
export  BBox_init, BBox_bound, BBox_expand, BBox_print, BBox_width
export  cardin, Dist, DistSq, VecPolygon2F
export  BBox_bottom_left, BBox_top_right


export  Polygon_prefix_lengths

export  Segment_nn_point
export  Segment_get_on, Segment_get_convex_coef
export  Segment_length
export  Segment_get_bisection_point

export  Polygon_length, Polygon_move_to_origin
export  Polygon_sample_uniformly, Polygon_push_smart, Polygon_spine

export  Polygon_read_plt_file


export  Polygon_simplify, Polygon_push, DistInfty
export  Polygon_simplify_radii
export  Polygon_simplify_ext
export  Polygon_translate!
export  Polygon_get_point_on
export  Polygon_as_matrix
export  Polygon_write_to_file

export  VecPnts_as_matrix

export  dist_seg_nn_point

export  convex_comb

end
