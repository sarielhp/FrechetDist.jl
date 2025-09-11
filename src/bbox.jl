## Axis aligned bounding box

module bbox

using StaticArrays
using Parameters
using LinearAlgebra
using DelimitedFiles
using Distributions
using Random

using ..cg.point
using ..cg.polygon

######################################################################
######################################################################
# Bounding box
"""
    BBox

Axis parallel bounding box.
"""

@with_kw_noshow mutable struct BBox{D,T}
    f_init::Bool = false
    mini::MVector{D,T} = zeros( T, D );
    maxi::MVector{D,T} = zeros( T, D );
end

BBox2F = BBox{2,Float64};

function  width( bb::BBox{D,T}, dim::Int64 = 1 ) where {D,T}
    return  bb.maxi[ dim ] - bb.mini[ dim ];
end
function  height( bb::BBox{D,T}, dim::Int64 = 2 ) where {D,T}
    return  bb.maxi[ dim ] - bb.mini[ dim ];
end


function Base.show(io::IO, bb::BBox{D,T}) where {D,T}
    for  i in 1:D
        print( io, " [[" );
        print( io, bb.mini[ i ] );
        print( io, ".." )
        print( io, bb.maxi[ i ] );
        print( io, "]] " );
    end
end

#=
function  print( bb::BBox{D,T}) where {D,T}
    for  i in 1:D
        Base.print( " [[" );
        Base.print( bb.mini[ i ] );
        Base.print( ".." )
        Base.print(  bb.maxi[ i ] );
        Base.print(  "]] " );
    end
end
=#

function  get_dim( bb, d )
    return  ( bb.mini[ d ], bb.maxi[ d ] )
end

function  middle( bb::BBox{D,T}, d ) where  {D,T}
    return   (bb.mini[ d ] + bb.maxi[ d ]) / 2.0;
end
function  d_min( bb::BBox{D,T}, d ) where  {D,T}
    return  bb.mini[ d ];
end
function  d_max( bb::BBox{D,T}, d ) where  {D,T}
    return  bb.maxi[ d ];
end

function  init( bb::BBox{D,T}, p, q ) where  {D,T}
    bb.f_init = true;
    for  i in 1:D
        bb.mini[ i ] = min( p[i ], q[ i ] );
        bb.maxi[ i ] = max( p[i ], q[ i ] );
    end
end

function  BBox2F_init( p::Point2F, q::Point2F )
    bb = BBox{2, Float64}()
    init( bb, p, q );
    return  bb
end

function  diam( bb::BBox{D,T} )  where  {D,T}
    ( ! bb.f_init )  &&  return  zero( T );

    sum = zero( T );
    for  i in 1:D
        sum += (bb.maxi[ i ] - bb.mini[ i ])^2;
    end

    return  sqrt( sum );
end

function  dist_in_dim( b::BBox{D,T}, c::BBox{D,T}, i::Int64
                             )  where  {D,T}

    ( b.mini[ i ] > c.maxi[ i ] )  &&  return  ( b.mini[ i ] - c.maxi[ i ] );
    ( b.maxi[ i ] < c.mini[ i ] )  &&  return  ( c.mini[ i ] - b.maxi[ i ] );

    return  zero( T );
end

"""
    extent_in_dim

    The maximum distance between two points in this dimension, that
belongs to the two bounding boxes.
"""
function  max_dist_in_dim( b::BBox{D,T}, c::BBox{D,T}, i::Int64
                               )  where  {D,T}

    ( b.mini[ i ] >= c.maxi[ i ] )  &&  return  ( b.maxi[ i ] - c.mini[ i ] );
    ( b.maxi[ i ] <= c.mini[ i ] )  &&  return  ( c.maxi[ i ] - b.mini[ i ] );

    # Boxes must intersect....
    return  max( abs( b.mini[ i ] - c.mini[ i ] ),
                 abs( b.maxi[ i ] - c.maxi[ i ] ),
                 abs( b.mini[ i ] - c.maxi[ i ] ),
                 abs( b.maxi[ i ] - c.mini[ i ] ) );
end

function  dist( b::BBox{D,T}, c::BBox{D,T} )  where  {D,T}
    ( ! b.f_init )  &&  return  zero( T );
    ( ! c.f_init )  &&  return  zero( T );

    sum = zero( T );
    for  i in 1:D
        sum += ( dist_in_dim( b, c, i ) )^2;
    end

    return  sqrt( sum );
end


function  max_dist( b::BBox{D,T}, c::BBox{D,T} )  where  {D,T}
    ( ! b.f_init )  &&  return  zero( T );

    sum = zero( T );
    for  i in 1:D
        sum += ( max_dist_in_dim( b, c, i ) )^2;
    end

    return  sqrt( sum );
end


function  bottom_left( bb::BBox{D,T} ) where  {D,T}
    return  Point{D,T}( bb.mini );
end

function  top_right( bb::BBox{D,T} ) where  {D,T}
    return  Point{D,T}( bb.maxi );
end

function  expand!( bb::BBox{D,T}, factor ) where  {D,T}
    mid::MVector{D,T} = (bb.maxi + bb.mini) / 2.0;
    diff::MVector{D,T} = (bb.maxi - bb.mini) * (factor/2.0);

#    margin::MVector{D,T} = (bb.maxi - bb.mini) * factor;
    bb.mini = mid - diff;
    bb.maxi = mid + diff;
end

function  expand( bb_in::BBox{D,T}, factor ) where  {D,T}
    bb = deepcopy( bb_in );
    expand!( bb, factor );
    return  bb;
end

function  expand_add!( bb::BBox{D,T}, _add ) where  {D,T}
    add::MVector{D,T} = fill( _add, D )
    mini::MVector{D,T} = bb.mini - add;
    maxi::MVector{D,T} = bb.maxi + add;

    bb.mini = mini;
    bb.maxi = maxi;
end

function  expand_add( _bb::BBox{D,T}, _add ) where  {D,T}
    bb = deepcopy( _bb );
    expand_add!( bb, _add );
    return  bb;
end

function Base.:+(b::BBox{D,T}, v::T ) where {D,T}
    return  expand_add( b, v );
end

function Base.in( p::Point{D,T}, bb::BBox{D,T} ) where {D,T}
    for  i ∈ 1:D
        if  ( ( bb.mini[ i ] > p[ i ] )  ||  ( bb.maxi[ i ] < p[ i ] ) )
            return  false
        end
    end

    return  true;
end

function  bound(  bb::BBox{D,T}, pnt::Point{D,T} )  where  {D,T}
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

function  bound(  bb::BBox{D,T}, P )  where  {D,T}
    for  p::Point{D,T}  in  P
        bound( bb, p )
    end
end

function  bound(  bb::BBox{D,T}, P::Polygon{D,T} )  where  {D,T}

    nv = cardin( P );
    for  i in 1:nv
        bound( bb, P[ i ] )
#        println( "i:", i, "   " );
#        print( bb );
#        println( "\n" );
    end
end

"""
    init

"""
function  init( P::Polygon{D,T},
                     range::UnitRange{Int64} ) where  {D,T}
    bb = BBox{D,T}()
    for  i in range
        bound( bb, P[ i ] );
    end

    return  bb;
end

function  init( P::Polygon{D,T} ) where  {D,T}
    bb = BBox{D,T}()
    for  p ∈ P
        bound( bb, p );
    end

    return  bb;
end

function  bound(  bb::BBox{D,T},
                       list::Vector{Polygon{D,T}} )  where  {D,T}
    for  x in list
        bound( bb, x )
    end
end


#####################################################################

function  top_left( bb::BBox2F )
    return  Point2F( bb.mini[ 1 ], bb.maxi[ 2 ] );
end
function  bottom_right( bb::BBox2F )
    return  Point2F( bb.maxi[ 1 ], bb.mini[ 2 ] );
end


#
export  expand, expand!
export  init,  bound, width
export  bottom_left, top_right
export  top_left, bottom_right
export  d_min, d_max, middle;
export  diam, dist;
export  dist, max_dist
export  width, height
export  BBox2F_init
export  BBox2F;
export  BBox;

end # // End module bbox

