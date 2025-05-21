module polygon

#using PrecompileTools
using Parameters
using StaticArrays
using LinearAlgebra
using DelimitedFiles

#include( "point.jl" );
using ..cg.point
#using point;
using ..cg.segment;

#import Point{D,T} = point.Point{D,T}

#include( "Polygon.jl" );

#using .Polygon;

### Conditional flag to determine whether a polygon is an array of
### points, or a struct containing a field which is an array of
### points. Seems to have no impact on performance. Arguably, having
### it as a distinct struct is "safer" type-wise.
const  POLYGON_AS_DIRECT_ARRAY = false;


###############################################33
###############################################33
### Polygon type
###############################################33

@static if  POLYGON_AS_DIRECT_ARRAY
    Polygon{D,T} = Vector{point.Point{D,T}};
else
    struct Polygon{D,T}
        pnts::Vector{point.Point{D,T}};
    end
end


@static if  POLYGON_AS_DIRECT_ARRAY
    @inline function Points( P::Polygon{D,T} ) where {D,T}
        return  P;
    end
else

    @inline function Points( P::Polygon{D,T} ) where {D,T}
        return  P.pnts;
    end
    @inline function Polygon{D,T}()  where {D,T}
        vec = Vector{Point{D,T}}();
        return  Polygon{D,T}( vec );
    end
end






function  Polygon_fill( P::Polygon{D,T}, f, range ) where {D,T}
    for  r ∈ range
        p = f( r );
        push!( P, Point{D,T}(p...) );
   end

    return  P;
end


function   Polygon_diff( src::Polygon{D,T}, sub::Polygon{D,T}
                  )::Polygon{D,T} where {D,T}
    PSX = setdiff( Points(src), Points( sub ) );
    PSA = Polygon_from_array( PSX );
    return  PSA;
end


function Polygon_from_array( arr::Vector{point.Point{D,T}} )  where {D,T}
    P = Polygon{D,T}();
    for  p in arr
        push!( P, p )
    end

    return  P;
end


function  Polygon_translate!( P::Polygon{D,T}, v::Point{D,T} ) where {D,T}
    for  i in 1:length(Points(P))
        # POINT p.x = p.x - v.x;
        P[i] = sub( P[i], v );
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
        push_smart!( Q, P[ t ] );
    end
    push_smart!( Q, (P[ i ] + P[ i+1 ] ) / 2.0 );
    for t in i+1:cardin(P)
        push_smart!( Q, P[ t ] );
    end
    return  Q;
end


function Base.ndims(::Type{Polygon{D, T}})::Int64 where {D,T}
    return  1;
end

function Base.show(io::IO, poly::Polygon{D,T}) where {D,T}
    f_iter::Bool = false;
    for p in polygon.Points( poly)
        if  f_iter
            print( io, "-" );
        else
            f_iter = true;
        end
        show(io, p );
    end
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
    m = zeros( T, D, length( Points(P) ) );
    for  i in 1:length( Points(P) )
        m[:,i] = P[ i ];
    end
    return  m;
end

"""
    write_to_file

    Writes a plt file of the polygon.
"""
function  write_to_file( P::Polygon{D,T}, filename ) where {D,T}
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

@static if  ! POLYGON_AS_DIRECT_ARRAY
    @inline function Base.last( poly::Polygon{D,T} ) where {D,T}
        #    println( "last??? " );
        return  last( Points( poly ) );
    end

    @inline function Base.first( poly::Polygon{D,T} ) where {D,T}
        #    println( "last??? " );
        return  first( Points( poly ) );
    end

    @inline function  Base.getindex(a::Polygon{D,T}, i::Int) where {D,T}
        return   a.pnts[ i ];
    end

    @inline function  Base.getindex(a::Polygon{D,T}, i) where {D,T}
        return   a.pnts[ i ];
    end

    @inline function  Base.setindex!(a::Polygon{D,T}, v, i::Int) where {D,T}
        a.pnts[ i ] = v;
    end
end

@static if  ! POLYGON_AS_DIRECT_ARRAY
    @inline function  Base.length( P::Polygon{D,T} ) where {D,T}
        return  length( Points(P) );
    end

    @inline function Base.eachindex( P::Polygon{D,T} ) where {D,T}
        return  eachindex( P.pnts );
    end

    @inline function  Base.iterate( P::Polygon{D,T} ) where {D,T}
        ( cardin(P) == 0 )  &&  return  nothing;
        return  (P[1], 1)
    end

    @inline function  Base.iterate( P::Polygon{D,T}, state ) where {D,T}
        ( state >= cardin( P) )  &&  return  nothing;

        return  (P[state+1], state+1)
    end
end


@static if  ! POLYGON_AS_DIRECT_ARRAY
    @inline function Base.push!( c::Polygon{D,T}, P::Polygon{D,T}) where {D,T}
        for p in P
            push!(c.pnts, p);
        end
    end

    @inline function Base.push!( c::Polygon{D,T}, p::Point{D,T}) where {D,T}
        push!(c.pnts, p );
    end

    @inline function Base.pop!( c::Polygon{D,T}) where {D,T}
        pop!(c.pnts);
    end
end

@inline function  Base.size( P::Polygon{D,T} ) where {D,T}
    return  Tuple{Int}( length( Points( P ) ) );
end

function  total_length( poly::Polygon{D,T} ) where {D,T}
    len::Float64 = 0.0;
    for i in firstindex(Points(poly)) : lastindex(Points(poly))-1
        len += Dist( poly[ i ], poly[ i + 1 ] );
    end

    return len;
end


function  Polygon_prefix_lengths( poly::Polygon{D,T}
)::Vector{Float64} where {D,T}
    v = Vector{Float64}();
    push!( v,  zero(Float64) );
    len::Float64 = 0.0;
    for i in firstindex(Points(poly)) : lastindex(Points(poly))-1
        p = poly[ i ];
        q = poly[ i + 1 ];
        len += Dist( p, q );
        push!( v, Float64( len ) );
    end
    return v;
end


@inline function  Polygon_edge_length( P::Polygon{D,T}, i::Int64 )  where {D,T}
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
        return first( Points(poly) );
    end
    if  t >= last( prefix_len )
        return  last( Points(poly) );
    end
    i = searchsortedfirst( prefix_len, t )
    #    println( "i= ", i );
    if  prefix_len[ i ] == t
        return  poly[ i ];
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

    return  convex_comb( poly[ i - 1 ], poly[ i ], x );
end


Polygon2I = Polygon{2,Int64};
Polygon2F = Polygon{2,Float64};

Polygon3F = Polygon{3,Float64};

@inline function  cardin( P::Polygon{D,T} )::Int64 where {D,T}
    return  length( Points( P ) );
end

"""
    Polygon_simplify_ext

Simplify the input curve P, where every vertex in distance larger than
r is added to the simplification. Maybe the cheapest and simplest
simplification algorithm. Works for all measures.

"""
function  Polygon_simplify_ext( P::Polygon{D,T}, r::T ) where {D,T}
#::Tuple{Polygon{D,T},Vector{Int64}}
    pout = Polygon{D,T}();
    pindices = Vector{Int64}();

    len = cardin( P );
    if  ( len == 0 )
        return pout;
    end
    if  ( push_smart!( pout, P[1] ) )
        push!( pindices, 1 );
    end

    curr = P[1];
    for i in 2:len
        if  ( Dist( P[i], curr ) > r )
            curr = P[ i ];
            if  ( push_smart!( pout, P[i] ) )
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

@inline function  push_smart!( pout::Polygon{D,T}, p::Point{D,T} ) where  {D,T}
    if  ( cardin( pout ) == 0 )  ||  ( DistSq( last( pout ), p ) > 0.0 )
        push!( pout, p );
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
        push_smart!( Q, P[ i ] );
        push_smart!( Q, p );
    end
    push_smart!( Q, P[ l ] );

    return  Q;
end


@inline function  Polygon_push( pout::Polygon{D,T}, p::Point{D,T} ) where  {D,T}
    push!( pout, p );
end



function Polygon_move_to_origin( P::Polygon{D,T} ) where {D,T}
    pout::Polygon{D,T} = Polygon{D,T}();
    p = P[ 1 ]
    for q in Points( P )
        push_smart!( pout, q - p );
    end
    return  pout;
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
    if  push_smart!( pout, P[1] )
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
            if  push_smart!( pout, P[i] )
                push!( pindices, i );
            end
        end
    end
    if  push_smart!( pout, P[ len ] )
        push!( pindices, len );
    end
    #println( "STRLEN: ", length( pindices ) )
    if  ( length( pindices ) <= 1 )
        Polygon_push( pout, P[ len ] )
        push!( pindices, len );
    end

    return  pout, pindices;
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
    push_smart!( new_P, first( Points(P) ) );

    sz = cardin( P );

    for  i  in  1:(sz-1)
        s = Segment{D,T}( P[ i ], P[ i + 1 ] );

        ell = Segment_length( s );

        ns::Int64 = floor( Int64, ell / delta );
        for  j  in 1:ns
            push_smart!( new_P, segment.at( s, j / ( ns + 1 ) ) );
        end
        push_smart!( new_P, P[ i + 1 ] );
    end

    return  new_P;
end


@inline function  spine( P::Polygon{D,T} )  where {D,T}
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

function  Polygon_random_gaussian( D,T,  n::Int64 )
    #x::T = zero(T);
    #x = x + x;

    P = Polygon{D,T}();
    for  i in 1:n
        p::Point{D,T} = Point_random_gaussian( D, T );
        push!( P, p  );
    end

    return P;
end


function  Polygon_random_sphere( D,T,  n::Int64 )
    #x::T = zero(T);
    #x = x + x;

    P = Polygon{D,T}();
    for  i in 1:n
        p::Point{D,T} = Point_random_gaussian( D, T );
        push!( P, normalize( p ) );
    end

    return P;
end


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
        push!( P, npoint( v[1], v[2] ) );
    end
    return P;
end


function  slice( P::Polygon{D,T}, d ) where {D,T}
    v = Vector{T}();
    for  p ∈ P  push!( v, p[ d ] ) end;
    return  v;
end


"""
    read_plt_orig_file

    Reads a .plt file into a polygon (2d floating point).
"""
function  read_plt_orig_file( filename, dchar = "," )
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
        push_smart!( P, npoint( x, y ) );
    end

    return  P
end

"""
    read_smf_file

Reads a .smf file into a polygon
"""
function  read_smf_file( filename )::Polygon3F
    P = Polygon3F();
    dchar = " "
    
    a = readdlm( filename );
    d_a::Int64 = size(a,1)

    for  r  in 1:size(a,1)
        prefix = a[r,1]
        #println( "line: [", prefix, "]" );
        if  prefix != "v"
            continue;
        end
        x,y,z = a[r, 2:4]
#        println( x, ", ", y, ", ",z );
        push_smart!( P, npoint( Float64(x), Float64(y), Float64(z) ) );
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


function  read_file( filename, dchar = "," )
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
        push_smart!( P, npoint( pieces[ 1 ], pieces[ 2 ] ) );
    end

    if  ( cardin( P ) <= 1 )
        Polygon_push( P, first( P ) );
    end
    return  P
end


function  read_txt_file( filename, dchar = "," )
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

function  times( P::Polygon{D,T} ) where {D,T}
    len = polygon.total_length( P );
    if  len == 0.0
        return  zeros( length( P ) );
    end
    push!( t, 0.0 );
    for  i  in  1:length(P) - 1
        delta = ( Dist( P[ i ], P[ i + 1] )
                  + Dist( Q[ i ], Q[ i + 1] ) ) / len;
        push!( t, last(t) + delta );
    end
    if last( t ) != 1.0
        pop!( t )
        push!( t, 1.0 );
    end
    return  t
end

function  at_time(
    P::Polygon{D,T},
    times::Vector{T},
    t::T
) where {D,T}
    ( t <= 0.0 ) &&  return first( P );
    ( t >= 1.0 ) &&  return last( P );

    pos = searchsortedfirst( times, t );
    @assert( 1 < pos <= length( P ) );
    prev = pos - 1;

    delta = ( t - times[ prev ] ) / (times[ pos ] - times[ prev ] );
    p = convex_comb( P[ prev ], P[ pos ], delta )

    return  p;
end


function  reverse( P::Polygon{D,T} ) where {D,T}
    Q = Polygon{D,T}();
    for  i ∈ length( P ):-1:1
        push!( Q, P[ i ] );
    end
    return  Q;
end

"""
    at( P, t)

Returns a point along P, parameterized uniformly on [0,1]. The
function at_time is faster, but requires preprocessing.
"""

function  at( P::Polygon{D,T}, t::T ) where {D,T}
    times = times( P );
    return  at_times( P, times, t );
end


function  (Base.:|>)( P::Polygon{D,T}, t::Tuple{Vararg{T,D}} )  where {D,T}
    push!( P, Point{D,T}( t... ) );
    return  P;
end


function  shift!( P::Polygon{D,T}, v::Point{D,T} ) where {D,T}
    for  i ∈ eachindex( P )
        P[i] = P[i] + v;
    end
    return  P;
end


function  mult!( P::Polygon{D,T}, scale::T ) where {D,T}
     for  i ∈ eachindex( P )
         P[i] = P[i] * scale;
     end
     return  P;
end

function  wiggle_seg( O::Polygon{D,T}, p, q, n, rate ) where {D,T}
    dir = q - p;
    vec = rate * dir;

    # Rotate by 90 degrees...
    #vec[ 2 ], vec[ 1 ] = vec[ 1 ], vec[ 2 ];
    if  D == 2
        vec = Point{D,T}( vec[ 2 ], vec[ 1 ] )
    else
        vec = Point{D,T}( vec[ 2 ], vec[ 1 ], vec[3:end]... )
    end

    delta = dir / ( n + 1 );
    curr = p;
    push_smart!( O, p );
    for  i  in 1:n
        curr = curr + delta;
        u = curr + vec;
        vec = -1.0 * vec;
        push_smart!( O, u );
    end
    push_smart!( O, q );
end


function  shortcut( P::Polygon{D,T}, i::Int64, j::Int64 ) where{D,T}
    Q = Polygon{D,T}();

    range = i+1:j-1;
    for  i ∈ 1:length( P )
        if  i ∈ range  continue  end;
        push!( Q, P[ i ] );
    end

    return  Q;
end


function  cut( P::Polygon{D,T}, r::UnitRange{Int} ) where{D,T}
    Q = Polygon{D,T}();

    for  i ∈ r
        push!( Q, P[ i ] );
    end

    return  Q;
end

function  append!( P::Polygon{D,T}, Q::Polygon{D,T} ) where{D,T}
    for  q ∈ Q
        push!( P, q );
    end
    
    return  P;
end

function  append_smart!( P::Polygon{D,T}, Q::Polygon{D,T} ) where{D,T}
    for  q ∈ Q
        push_smart!( P, q );
    end
    
    return  P;
end


"""
    wiggle

Generates a version of the original polygon, where each segment is
replaced by a "wiggle" with n vertices, and displacement of
wave_size. Intended for 2d, but might work in higher dimensions.

"""

function  wiggle( P::Polygon{D,T}, n::Int64, wave_size::Float64 ) where {D,T}
    O = Polygon{D,T}();
    ( total_length( P ) == 0.0 )  &&  return  O;

    for i ∈ 1:length(P) - 1
        wiggle_seg( O, P[i], P[i + 1], n, wave_size );
    end

    return  O;
end

####################n#######################################################
###########################################################################
###########################################################################


export  Polygon_move_to_origin
export  Polygon_sample_uniformly, push_smart!, spine

export  read_file
export  read_plt_orig_file
export  read_txt_file
export  read_smf_file
export  write_to_file

export  Polygon_prefix_lengths
export  Polygon_simplify, Polygon_push, DistInfty
export  Polygon_simplify_radii
export  Polygon_translate!
export  Polygon_get_point_on
export  Polygon_as_matrix

export  Polygon_shortcut
export  Polygon_random
export  Polygon_random_gaussian
export  Polygon_random_sphere

export  Polygon_convex_comb
export  Polygon_split_edges
export  Polygon_split_single_edge
export  Polygon_edge_length
export  Polygon_from_array
export  Polygon_sample
export  Polygon_diff
export  Polygon_fill

export  total_length;

export  Polygon;
export  Polygon2I, Polygon2F;
export  cardin
export  Polygon_simplify_ext
export  slice
export  reverse, Points

export  VecPnts_as_matrix

export  times, at_time, at;
export  shift!, mult!, ndims, cut, append!, append_smart!

end # // End module polygon
