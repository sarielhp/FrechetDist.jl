# Palette
#
# A palette of a polygon P[1..n] is an array A[1..n], where A[i] is the
# Frechet error introduced by the ith vertex. The array is filled in a
# binary search hierarchical way, so that all the vertices with error
# above a certain threshold can be extracted in time linear in their
# number, instead of linear time in the number of vertices of the
# curve. The construction and extraction are generic, and live in
# BasicCompGeometry.PntSeqHierarchy. Here we only supply the Frechet
# error of a subcurve, and the midpoint splitting rule.

const PH = BasicCompGeometry.PntSeqHierarchy

"""
    frechet_palette

The palette of `P`, as a vector with one entry per vertex.
"""
frechet_palette( P::PntSeq{N,T} ) where {N,T} =
    PH.palette( P, frechet_width_approx ).err::Vector{T}

"""
    frechet_palette_approx

Extract an approximation to P using the palette, and a width
parameter. The output polygon Q has the property the Frechet distance
betwwen P and Q is at most w (might be smaller).

The approximation provided is fast (i.e., running time proportional to
the nubmer of vertices of Q, but the approximation might be of low
quality.

"""
frechet_palette_approx( P::PntSeq{N,T}, plt::Vector{T}, w::T ) where {N,T} =
    PH.palette_extract( P, plt, w )

"""
    frechet_palette_level

Extract the approximate polygon by specifying a level of the
approximation.

"""
function frechet_palette_level( P::PntSeq{N,T}, plt::Vector{T},
                                level::Int64 ) where {N,T}
    pout, pout_is = PH.palette_level( P, plt, level );
    println( "|pout|: ", cardin( pout ) );
    return  pout, pout_is;
end
