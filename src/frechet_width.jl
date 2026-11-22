# Fast estimates of the Frechet distance between a subcurve and the segment joining
# its endpoints (the width), used by the simplifiers and for vertex offsets.

"""
    frechet_width_approx

 2-approximation to the Frechet distance between
   seg = P[first(rng)]-P[last(rng)] and he polygon
   P[rng]
 Here, rng is a range i:j

This function implements a greedy matching alnog the segment - you
move to the cloest point on seg to the current point, making sure
never moving back.

"""
function frechet_width_approx(P::PntSeq{N,T}, rng::UnitRange{Int64} = 0:0)::T where {N,T}
    check_coordinate_type(P, "frechet_width_approx")
    card = cardin(P)
    if (card <= 2)
        return 0
    end

    if (rng == 0:0)
        rng = 1:card
    end
    (isempty(rng) || (1 <= first(rng) && last(rng) <= card)) || throw(
        ArgumentError("frechet_width_approx: range $rng is outside the vertices 1:$card"),
    )

    if (length(rng) <= 2)
        return 0
    end

    return sqrt(width_leash_sq(P, rng, Inf))
end

# The squared greedy width of P[rng] against its spine (rng has at least 3 vertices).
# Stops early once sqrt of the running maximum exceeds `w`: the maximum only grows and
# sqrt is monotone, so the comparison with `w` is decided at that point.
function width_leash_sq(P::PntSeq{N,T}, rng::UnitRange{Int64}, w::Float64) where {N,T}
    s_p = P[first(rng)]
    s_q = P[last(rng)]

    seg_len_sq = fp_dist_sq(s_p, s_q)
    s_vec = s_q - s_p

    t::Float64 = 0
    curr::Point{N,T} = s_p
    leash::Float64 = 0
    for i in (first(rng)+1):(last(rng)-1)
        q, new_t = fp_nn_point_ext_ext(s_p, s_q, P[i], seg_len_sq, s_vec)
        if (new_t > t)
            t = new_t
            curr = q
        end
        leash = max(leash, fp_dist_sq(curr, P[i]))
        (sqrt(leash) > w) && return leash
    end
    return leash
end

"""
    frechet_width_exceeds( P, rng, w )

Whether `frechet_width_approx(P, rng) > w`, computed with the same arithmetic but
stopping as soon as the answer is known.
"""
function frechet_width_exceeds(
    P::PntSeq{N,T},
    rng::UnitRange{Int64},
    w::Float64,
) where {N,T}
    (cardin(P) <= 2 || length(rng) <= 2) && return zero(T) > w
    return sqrt(width_leash_sq(P, rng, w)) > w
end

eachindexButLast(x) = firstindex(x):(lastindex(x)-1)

#####################################################################
"""
    frechet_offsets( P, p_indices )

Given `P` and the increasing indices `p_indices` of a subcurve of `P`, returns, for
every edge `p_indices[i] p_indices[i+1]` of the subcurve, a fast upper estimate
(`frechet_width_approx`) of the Frechet distance between that edge and the portion
of `P` it replaces: one value per edge.

"""
function frechet_offsets(P::PntSeq{N,T}, p_indices::Vector{Int64}) where {N,T}
    check_coordinate_type(P, "frechet_offsets")
    check_indices(p_indices, cardin(P), "frechet_offsets")
    offs = Vector{Float64}()

    for i in eachindexButLast(p_indices)
        d = frechet_width_approx(P, p_indices[i]:p_indices[i+1])
        push!(offs, d)
    end
    return offs
end
