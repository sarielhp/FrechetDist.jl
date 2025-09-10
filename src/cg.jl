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
#include

include( "segment.jl" );
using .segment

include( "polygon.jl" );
using .polygon

include( "trans2d.jl" );
using .trans2d

include( "bbox.jl" );
using .bbox


include( "polygon_hausdorff.jl" );
using .polygon_hausdorff



using LinearAlgebra
using DelimitedFiles
#using .point

#using .Point;
#using .Polygon;
#using LoopVectorization
#using SIMD

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



distance(p1::Point{D,T}, p2::Point{D,T}) where {D,T} = norm( sub( p1, p2) )

function distance(y::Point{D,T}, l::Line{D,T}) where {D,T}
    p, u = l.p, l.u

    t = dot( (y - p), u ) / dot( u, u )
    x = Point(p + t*u)

    return Dist(x, y)
end


#######################################################################3


#####################################################################
# Predefined useful types...




VecPoint2I = Vector{Point{2,Int64}};
VecPolygon2F = Vector{Polygon{2,Float64}};



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

export   is_left_turn
export   is_right_turn




#export  iseg_nn_point
#export  iseg_nn_point_ext

export  iseg_iseg_dist


export  Points;
export  Centroid;

export  convex_comb

export  segs_match_price
export  cardin, VecPolygon2F


end
