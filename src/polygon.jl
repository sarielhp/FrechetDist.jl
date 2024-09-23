module polygon

#using PrecompileTools
using Parameters
using StaticArrays
using LinearAlgebra
using DelimitedFiles

#include( "point.jl" );
using point;
using segment;

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
    function Points( P::Polygon{D,T} ) where {D,T}
        return  P;
    end
else
    
    function Points( P::Polygon{D,T} ) where {D,T}
        return  P.pnts;
    end
    function Polygon{D,T}()  where {D,T}
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
    function Base.last( poly::Polygon{D,T} ) where {D,T}
        #    println( "last??? " );
        return  last( Points( poly ) );
    end

    function Base.first( poly::Polygon{D,T} ) where {D,T}
        #    println( "last??? " );
        return  first( Points( poly ) );
    end

    function  Base.getindex(a::Polygon{D,T}, i::Int) where {D,T}
        return   a.pnts[ i ];
    end

    function  Base.setindex!(a::Polygon{D,T}, v, i::Int) where {D,T}
        a.pnts[ i ] = v;
    end
end

@static if  ! POLYGON_AS_DIRECT_ARRAY
    function  Base.length( P::Polygon{D,T} ) where {D,T}
        return  length( Points(P) );
    end

    function Base.eachindex( P::Polygon{D,T} ) where {D,T}
        return  eachindex( P.pnts );
    end

    function  Base.iterate( P::Polygon{D,T} ) where {D,T}
        ( cardin(P) == 0 )  &&  return  nothing;
        return  (P[1], 1)
    end

    function  Base.iterate( P::Polygon{D,T}, state ) where {D,T}
        ( state >= cardin( P) )  &&  return  nothing;

        return  (P[state+1], state+1)
    end
end


@static if  ! POLYGON_AS_DIRECT_ARRAY
    function Base.push!( c::Polygon{D,T}, P::Polygon{D,T}) where {D,T}
        for p in P
            push!(c.pnts, p);
        end
    end

    function Base.push!( c::Polygon{D,T}, p::Point{D,T}) where {D,T}
        push!(c.pnts, p);
    end

    function Base.pop!( c::Polygon{D,T}) where {D,T}
        pop!(c.pnts);
    end
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

function  cardin( P::Polygon{D,T} )::Int64 where {D,T}
    return  length( Points( P ) );
end

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

function  push_smart!( pout::Polygon{D,T}, p::Point{D,T} ) where  {D,T}
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
        push_smart!( Q, P[ i ] );
        push_smart!( Q, p );
    end
    push_smart!( Q, P[ l ] );

    return  Q;
end


function  Polygon_push( pout::Polygon{D,T}, p::Point{D,T} ) where  {D,T}
    push!( pout, deepcopy( p ) );
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
            push_smart!( new_P,
                                Segment_get_on( s, j / ( ns + 1 ) ) );
        end
        push_smart!( new_P, P[ i + 1 ] );
    end

    return  new_P;
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




export  Polygon_move_to_origin
export  Polygon_sample_uniformly, push_smart!, Polygon_spine

export  read_file
export  read_plt_orig_file
export  read_txt_file
export  write_to_file

export  Polygon_prefix_lengths
export  Polygon_simplify, Polygon_push, DistInfty
export  Polygon_simplify_radii
export  Polygon_translate!
export  Polygon_get_point_on
export  Polygon_as_matrix
export  Polygon_random
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

export  VecPnts_as_matrix


end # // End module polygon