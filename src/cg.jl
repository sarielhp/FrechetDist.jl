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


include( "polygon_hausdorff.jl" );
using .polygon_hausdorff



using LinearAlgebra
using DelimitedFiles
#using .point

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
    ( ! b.f_init )  &&  return  zero( T );
    ( ! c.f_init )  &&  return  zero( T );

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

function  BBox_expand!( bb::BBox{D,T}, factor ) where  {D,T}
    mid::MVector{D,T} = (bb.maxi + bb.mini) / 2.0;
    diff::MVector{D,T} = (bb.maxi - bb.mini) * (factor/2.0);

#    margin::MVector{D,T} = (bb.maxi - bb.mini) * factor;
    bb.mini = mid - diff;
    bb.maxi = mid + diff;
end

function  BBox_expand( bb_in::BBox{D,T}, factor ) where  {D,T}
    bb = deepcopy( bb_in );
    BBox_expand!( bb, factor );
    return  bb;
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

function  BBox_bound(  bb::BBox{D,T}, P )  where  {D,T}
    for  p::Point{D,T}  in  P
        BBox_bound( bb, p )
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

function  BBox_init( P::Polygon{D,T} ) where  {D,T}
    bb = BBox{D,T}()
    for  p ∈ P
        BBox_bound( bb, p );
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




VecPoint2I = Vector{Point{2,Int64}};
VecPolygon2F = Vector{Polygon{2,Float64}};
BBox2F = BBox{2,Float64};

#####################################################################


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


#
export  BBox_expand, BBox_expand!
export  BBox_init, BBox_bound, BBox_print, BBox_width
export  cardin, VecPolygon2F
export  BBox_bottom_left, BBox_top_right
export  BBox_min, BBox_max, BBox_middle;
export  BBox_diam, BBox_dist;
export  BBox_dist, BBox_max_dist

export   BBox_width, BBox_height;


#export  iseg_nn_point
#export  iseg_nn_point_ext

export  iseg_iseg_dist


export  Points;
export  Centroid;

export  convex_comb

export  segs_match_price

end
