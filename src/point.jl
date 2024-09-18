module point

using StaticArrays
using Parameters
using LinearAlgebra
using DelimitedFiles


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
Point2I = Point{2,Int64};
Point2F = Point{2,Float64};

#    x::MVector{D,T} = zeros( T, D );
#    return  Point{D,T}( x );
#=
function Point{D,T}() where  {D,T}
    return  Point{D,T}( zeros( T, D ) );
end
=#

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
    @inbounds for i in 1:D
        @fastmath sum +=  p[i]^2;
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
#    return  sum( p .* q );
    s = zero( T );
    @inbounds for i in 1:D
        s += p[ i ] * q[ i ];
    end
    return  s;

    #return  LinearAlgebra.dot( p, q );
end

#=
# Define standard lexicographical ordering on points
function Base.isless( p::Point{D,T}, q::Point{D,T} ) where {D,T}
    @inbounds for i in 1:D
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
=#

function  DistSq(p::Point{D,T}, q::Point{D,T})::T where {D,T}
    sum = zero(T);
    @inbounds for  i in 1:D
        @fastmath sum += ( p[i] - q[i] )^2
    end

    return  sum  #norm( p1.x - p2.x );
end

function  Dist(p::Point{D,T}, q::Point{D,T}) where {D,T}
    sum = 0.0;
    @inbounds for  i in 1:D
        @fastmath sum +=  ( p[i] - q[i] )^2
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
    o = Point{D,T}(undef )

    #@inbounds
    @inbounds for  i in 1:D
        #o[ i ] = p[ i ] * (1.0 - t)  + q[ i] * t;
        @fastmath o[ i ] = p[ i ]  + ( q[ i] - p[ i ] ) * t;
#        o[ i ] = p[ i ] * s + q[ i] * t;
    end

    return  o;
#    return  add( mult( p, 1.0-t), mult( q, t ) );
end


function  is_left_turn( p::Point2F, q::Point2F, r::Point2F )
    q_x = q[1] - p[1];
    q_y = q[2] - p[2];
    r_x = r[1] - p[1];
    r_y = r[2] - p[2];

    return   ( q_x * r_y - q_y * r_x ) >  0.0;
end
function  is_right_turn( p::Point2F, q::Point2F, r::Point2F )
    q_x = q[1] - p[1];
    q_y = q[2] - p[2];
    r_x = r[1] - p[1];
    r_y = r[2] - p[2];

    return   ( q_x * r_y - q_y * r_x ) <  0.0;
end

"""
    npoint( args... )

A flexible constructor for a point specified by the arguments. Thus
npoint( 2.0, 3.0, 4.0 ) defined the 3d point (2.0, 3.0, 4.0). Or
similarly, point( 2.0, 1.0 ) would create a 2d point.

"""
function npoint( args...)
    D=length(args);
    T=typeof( first( args ) )
    x::MVector{D,T} = MVector{D,T}(undef);

    for i in eachindex( x )
        x[ i ] = args[ i ];
    end

    p::Point{D,T} = Point{D,T}( x );

    return  p;
end


function  Point_random( D,T )
    x::MVector{D,T} = MVector{D,T}( undef );

    @inbounds for i in eachindex( x )
        x[ i ] = T( rand() );
    end

    p::Point{D,T} = Point{D,T}( x );

    return  p;
end



export Point
export Point2F

export DistSq
export Dist
export Point_random, sub, convex_comb, dot

export  add, sub, mult, norm

export npoint

end # end to module Point
