module trans2d

using Parameters
using StaticArrays
using LinearAlgebra
using DelimitedFiles

include("point.jl")
using .point;

#import Base:


function  translation( a::T, b::T ) where {T}
    return  [ one(T) zero(T) a; zero(T) one(T) b; zero(T) zero(T) one(T) ];
end

function  translation( p::Point{2,T} ) where {T}
    return  translation( p[1], p[2] );
end

function  rotation( x::T ) where {T}
    return [   cos(x)  sin(x)   zero(T);
              -sin(x)  cos(x)   zero(T);
              zero(T)  zero(T)  one(T) ];
end


function  mult( m::Matrix{T}, P::Point{2,T} )::Point{2,T}  where {T}
    res = m * [ P[1], P[2], 1.0 ]
    return  Point{2,T}( res[ 1 ], res[ 2 ] );
end


export  translation, rotation, mult

end
