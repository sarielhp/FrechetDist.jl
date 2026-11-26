# Palette
#
# A palette of a polygon P[1..n] is an array A[1..n], where A[i] is the
# Frechet error introduced by the ith vertex. The array is filled in a
# hierarchical binary-search way, so that all the vertices with error
# above a certain threshold can be extracted in time linear in their
# number, instead of time linear in the number of vertices of the
# curve. The construction and extraction are generic, and live in
# BasicCompGeometry.PntSeqHierarchy. Here we only supply the Frechet
# error of a subcurve, and the midpoint splitting rule.

const PH = BasicCompGeometry.PntSeqHierarchy

"""
    frechet_palette

The palette of `P`, as a vector with one entry per vertex.
"""
function frechet_palette(P::GPntSeq{N,T,PT}) where {N,T,PT<:AbstractVector{T}}
    check_coordinate_type(P, "frechet_palette")
    return PH.palette(P, frechet_width_approx).err::Vector{T}
end

"""
    frechet_palette_approx

Extract an approximation to P using the palette, and a width
parameter. The output polygon Q has the property that the Frechet distance
between P and Q is at most w (it might be smaller).

The approximation provided is fast (i.e., running time proportional to
the number of vertices of Q), but the approximation might be of low
quality.

"""
function frechet_palette_approx(
    P::GPntSeq{N,T,PT},
    plt::Vector{T},
    w::Real,
) where {N,T,PT<:AbstractVector{T}}
    check_palette(P, plt, "frechet_palette_approx")
    return PH.palette_extract(P, plt, w)
end

function check_palette(P::PntSeq, plt::Vector, who::AbstractString)
    check_coordinate_type(P, who)
    length(plt) == cardin(P) || throw(
        ArgumentError(
            "$who: the palette has $(length(plt)) entries; the curve has $(cardin(P)) vertices",
        ),
    )
    return nothing
end

"""
    frechet_palette_level

Extract the approximate polygon by specifying a level of the
approximation.

"""
function frechet_palette_level(
    P::GPntSeq{N,T,PT},
    plt::Vector{T},
    level::Int64,
) where {N,T,PT<:AbstractVector{T}}
    check_palette(P, plt, "frechet_palette_level")
    pout, pout_is = PH.palette_level(P, plt, level)
    return pout, pout_is
end
