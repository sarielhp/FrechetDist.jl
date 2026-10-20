###########################################################################
###########################################################################
###########################################################################
# Greedy simplification


function find_frechet_prefix_inner(
    P::PntSeq{D,T},
    strt::Int64,
    i::Int64,
    j::Int64,
    w::T,
) where {D,T}
    if (i >= j) || ((strt + 1) == j)
        return j # the width is zero, nothing to do.
    end
    r = frechet_width_approx(P, strt:j)
    if (r <= w)
        return j
    end
    if (i+1) == j
        return i
    end
    j = j - 1
    mid = (i + j) >> 1
    r_m = frechet_width_approx(P, strt:mid)
    if (r_m > w)
        return find_frechet_prefix_inner(P, strt, i, mid - 1, w)
    end

    # No reason to waste time?
    if (j - mid) <= ((j - strt) ÷ 10)
        #println( "FLOGI" );
        return mid
    end

    return find_frechet_prefix_inner(P, strt, mid, j, w)
end

function find_frechet_prefix(P::PntSeq{D,T}, i::Int64, j::Int64, w::T) where {D,T}
    r = frechet_width_approx(P, i:j)
    if (r <= w)
        return j
    end
    return find_frechet_prefix_inner(P, i, i, j, w)
end

function frechet_simplify_to_width(P::PntSeq{D,T}, w::T) where {D,T}
    check_curve(P, "frechet_simplify_to_width", "the curve"; min_vertices = 0)
    check_nonnegative(w, "frechet_simplify_to_width", "w")
    pout = PntSeq{D,T}()
    pindices = Vector{Int64}()

    len = cardin(P)
    if (len == 0)
        return pout, pindices
    end
    if (push_distinct!(pout, P[1]))
        push!(pindices, 1)
    end

    curr_ind = 1
    next_ind::Int64 = 1
    card = cardin(P)
    while true
        next_ind = find_frechet_prefix(P, curr_ind, card, w)
        #println( "next_ind: ", next_ind );
        @assert(next_ind > curr_ind)
        push!(pindices, next_ind)
        push!(pout, P[next_ind])
        if next_ind == card
            return pout, pindices
        end
        curr_ind = next_ind
    end
end


"""
    frechet_simplify_to_cardin( P, sz )

Simplifies `P` by refining the allowed Frechet error until the simplification has
at least `sz + 2` vertices, or at least half the vertices of `P`, whichever comes
first. The vertices kept are a subset of those of `P`, in order. Returns the
simplification and the indices of the kept vertices in `P`. Note that the number of
vertices returned can be considerably larger than `sz`.

When the simplification found is exact (Frechet distance zero, for example for a
straight track), no curve with more vertices exists at that accuracy, and it is
returned even if it has fewer than `sz + 2` vertices.
"""
function frechet_simplify_to_cardin(P::PntSeq{D,T}, sz::Int64) where {D,T}
    check_curve(P, "frechet_simplify_to_cardin", "the curve")
    pout = PntSeq{D,T}()
    pindices = Vector{Int64}()

    @assert(cardin(P) >= 2)
    if sz <= 2
        push!(pindices, 1, cardin(P))
        return spine(P), pindices
    end
    w = frechet_width_approx(P)
    ratio = 0.8

    count::Int64 = 1
    while true
        count = count + 1
        w = w * ratio
        if count > 5
            w = w / count
        end
        Q, q_indices = Polygon_simplify_ext(P, w)
        m = frechet_c_mono_approx_subcurve(P, Q, q_indices)[1]
        R, R_indices = frechet_simplify_to_width(P, m.leash)
        if ((cardin(R) >= sz + 2) || (2*cardin(R) >= cardin(P)))
            return R, R_indices
        end
        m = frechet_c_mono_approx_subcurve(P, R, R_indices)[1]
        # R is already exact: a smaller width cannot produce a larger simplification
        # (without this the loop below never ends, since w stays 0).
        (m.leash == 0.0) && return R, R_indices
        w = m.leash
    end
end


###########################################################################



function exp_search_width_prefix(P::PntSeq{N,T}, start::Int64, w::Float64) where {N,T}
    len = cardin(P)
    hi::Int64 = min(start + 2, len)

    (hi >= len) && return hi
    while hi < len
        r = frechet_width_approx(P, start:hi)
        (r > w) && return min(hi + 5, len)
        hi = start + 2 * (hi - start)
    end

    return len
end


function frechet_simplify_w_exp(P::PntSeq{D,T}, w::T) where {D,T}
    check_curve(P, "frechet_simplify_w_exp", "the curve"; min_vertices = 0)
    check_nonnegative(w, "frechet_simplify_w_exp", "w")
    pout = PntSeq{D,T}()
    pindices = Vector{Int64}()

    card = cardin(P)
    if (card == 0)
        return pout, pindices
    end
    if (push_distinct!(pout, P[1]))
        push!(pindices, 1)
    end

    curr_ind = 1
    next_ind::Int64 = 1
    while true
        hi = exp_search_width_prefix(P, curr_ind, w)
        next_ind = find_frechet_prefix(P, curr_ind, hi, w)
        @assert(next_ind > curr_ind)
        push!(pindices, next_ind)
        push!(pout, P[next_ind])
        if next_ind == card
            return pout, pindices
        end
        curr_ind = next_ind
    end
end


#
# End of file
##########################################################################
