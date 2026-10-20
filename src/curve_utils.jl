# Curve helpers specific to the Frechet distance code. Everything that is
# generic geometry (points, segments, point sequences, bounding boxes)
# lives in BasicCompGeometry.

"""
    Polygon_simplify_ext

Simplify the input curve P, where every vertex in distance larger than
r is added to the simplification. Maybe the cheapest and simplest
simplification algorithm. Works for all measures. Returns the
simplification and the indices of the vertices kept (both empty for an empty curve).
"""
function Polygon_simplify_ext(P::PntSeq{D,T}, r::T) where {D,T}
    pout = PntSeq{D,T}()
    pindices = Vector{Int64}()

    len = cardin(P)
    if (len == 0)
        return pout, pindices
    end
    if (push_distinct!(pout, P[1]))
        push!(pindices, 1)
    end

    curr = P[1]
    for i = 2:len
        if (fp_dist(P[i], curr) > r)
            curr = P[i]
            if (push_distinct!(pout, P[i]))
                push!(pindices, i)
            end
        end
    end

    # We push the last vertex in event if it is a repetition, for
    # example if the curve is close, for example...
    push!(pout, P[len])
    push!(pindices, len)

    return pout, pindices
end

####################################################################
# We specify for each vertex how much error it is willing to accept
# What we do here is pretty Conservative and silly. Should be enough
# as a first try.
####################################################################
function Polygon_simplify_radii(P::PntSeq{D,T}, r::Vector{T}) where {D,T}
    @assert(cardin(P) == length(r))

    pindices = Vector{Int64}()
    pout = PntSeq{D,T}()

    len = cardin(P)
    if (len == 0)
        return pout, pindices
    end
    if push_distinct!(pout, P[1])
        push!(pindices, 1)
    end

    curr = P[1]
    curr_r = r[1]
    for i = 2:len
        curr_r = min(curr_r, r[i])
        if (fp_dist(P[i], curr) > curr_r)
            curr = P[i]
            if (i < len)
                curr_r = r[i+1]
            end
            if push_distinct!(pout, P[i])
                push!(pindices, i)
            end
        end
    end
    if push_distinct!(pout, P[len])
        push!(pindices, len)
    end
    if (length(pindices) <= 1)
        push!(pout, P[len])
        push!(pindices, len)
    end

    return pout, pindices
end

############################################################################
# Returns a polygon where the points are sampled uniformly along the polygon
# n is roughly the number of vertices one gets
############################################################################
function Polygon_sample_uniformly(P::PntSeq{D,T}, n::Int64) where {D,T}
    len = fp_geom_length(P)
    delta = len / (n-1)

    new_P = PntSeq{D,T}()
    push_distinct!(new_P, first(P))

    for i = 1:(cardin(P)-1)
        s = Segment{D,T}(P[i], P[i+1])

        ns::Int64 = floor(Int64, fp_geom_length(s) / delta)
        for j = 1:ns
            push_distinct!(new_P, fp_at(s, j / (ns + 1)))
        end
        push_distinct!(new_P, P[i+1])
    end

    return new_P
end

"""
    remove_consecutive_duplicates( P )

The curve `P` without consecutive repeated vertices. Returns `P` itself when there
are none, or when fewer than two distinct vertices would remain (a curve needs two
vertices).
"""
function remove_consecutive_duplicates(P::PntSeq{D,T}) where {D,T}
    Q = PntSeq{D,T}()
    for p in P
        push_distinct!(Q, p)
    end
    (cardin(Q) < 2 || cardin(Q) == cardin(P)) && return P
    return Q
end

"""
    spine

The curve made out of the first and last vertex of `P`.
"""
function spine(P::PntSeq{D,T}) where {D,T}
    pout = PntSeq{D,T}()
    if (cardin(P) <= 0)
        return pout
    end
    push!(pout, P[1])
    push!(pout, last(P))

    return pout
end

"""
    read_file

Reads a curve from a text file with one point per line (see
`BasicCompGeometry.ReadWrite.read_plt`). A curve with at most one point is
completed to a two point curve by repeating that point.
"""
function read_file(filename, dchar = ",")
    P = read_plt(filename)
    if (cardin(P) <= 1)
        push!(P, first(P))
    end
    return P
end
