# Geometric primitives with the floating-point handling the Frechet code
# relies on. They differ from the BasicCompGeometry versions in rounding
# (BCG uses a different formulation) and in snapping: results that are
# within 1e-6 of an endpoint are returned as the endpoint itself, and
# degenerate segments are treated explicitly. The Frechet code depends on
# these exact values, so they stay local rather than using BCG's.

@inline function fp_dot(p::Point{D,T}, q::Point{D,T}) where {D,T}
    s = zero(T)
    @inbounds for i = 1:D
        s += p[i] * q[i]
    end
    return s
end

@inline function fp_norm(p::Point{D,T}) where {D,T}
    sum = zero(T)
    @inbounds for i = 1:D
        @fastmath sum += p[i]^2 # julia-audit: allow-fastmath (bit-identical to the original numerics)
    end
    return sqrt(sum)
end

@inline function fp_dist_sq(p::Point{D,T}, q::Point{D,T})::T where {D,T}
    sum = zero(T)
    @inbounds for i = 1:D
        @fastmath sum += (p[i] - q[i])^2 # julia-audit: allow-fastmath (bit-identical to the original numerics)
    end
    return sum
end

@inline function fp_dist(p::Point{D,T}, q::Point{D,T}) where {D,T}
    sum = 0.0
    @inbounds for i = 1:D
        @fastmath sum += (p[i] - q[i])^2 # julia-audit: allow-fastmath (bit-identical to the original numerics)
    end
    return sqrt(sum)
end

"""
    fp_convex_comb( p, q, t )

`(1-t) p + t q`, returning `p` and `q` exactly at `t == 0` and `t == 1`. There is no
snapping band: a fixed band in `t` moves points by a distance proportional to the edge
length, which distorts long edges.
"""
@inline function fp_convex_comb(p::Point{D,T}, q::Point{D,T}, t::Float64) where {D,T}
    (t == 0.0) && return p
    (t == 1.0) && return q

    o = MVector{D,T}(undef)
    @inbounds for i = 1:D
        @fastmath o[i] = p[i] + (q[i] - p[i]) * t # julia-audit: allow-fastmath (bit-identical to the original numerics)
        @assert(! isnan(o[i]))
    end

    return Point{D,T}(o)
end

@inline fp_at(seg::Segment{D,T}, t::Float64) where {D,T} = fp_convex_comb(seg.p, seg.q, t);

"""
    fp_nn_point_ext_ext

Nearest point on the segment `s_p`-`s_q` to `qr`, together with its convex
combination coefficient. The caller must precompute

    seg_len_sq = fp_dist_sq( s_p, s_q )
    s_vec      = s_q - s_p
"""
@inline function fp_nn_point_ext_ext(
    s_p::Point{D,T},
    s_q::Point{D,T},
    qr::Point{D,T},
    seg_len_sq::T,
    s_vec::Point{D,T},
) where {D,T}
    # v(t) = p + t*(q-p),  D(t) = ||p-qr||^2 + 2t <p-qr,q-p> + t^2 ||q-p||^2
    # The minimum is achieved at  t^* = -b/(2a),  a = ||q-p||^2,
    # b = 2 <p-qr,q-p>.
    a = seg_len_sq
    if (a == 0.0)
        return s_p, zero(T)
    end
    b = 2.0 * fp_dot(s_p - qr, s_vec)
    t::T = -b / (2.0 * a)

    if (t <= 0.0)
        return s_p, zero(T)
    elseif (t >= 1.0)
        return s_q, one(T)
    end

    return fp_convex_comb(s_p, s_q, t), t
end

"""
    fp_nn_point_ext

As `fp_nn_point_ext_ext`, without the precomputed parameters.
"""
@inline function fp_nn_point_ext(
    s_p::Point{D,T},
    s_q::Point{D,T},
    qr::Point{D,T},
) where {D,T}
    return fp_nn_point_ext_ext(s_p, s_q, qr, fp_dist_sq(s_p, s_q), s_q - s_p)
end

"""
    fp_nn_point

Nearest point on the segment `s_p`-`s_q` to the query point `qr`.
"""
@inline function fp_nn_point(s_p::Point{D,T}, s_q::Point{D,T}, qr::Point{D,T}) where {D,T}
    a = fp_dist_sq(s_p, s_q)
    if (a == 0.0)
        return s_p
    end
    b = 2.0 * fp_dot(s_p - qr, s_q - s_p)
    t::T = -b / (2.0 * a)

    if (t <= 0.0)
        return s_p
    end
    if (t >= 1.0)
        return s_q
    end

    return fp_convex_comb(s_p, s_q, t)
end
@inline fp_nn_point(s::Segment{D,T}, qr::Point{D,T}) where {D,T} = fp_nn_point(s.p, s.q, qr);

"""
    fp_convex_coef( seg, p )

Coefficient `d(seg.p, p) / |seg|`, or 0 for a degenerate segment.
"""
@inline function fp_convex_coef(s::Segment{D,T}, qr::Point{D,T}) where {D,T}
    len = fp_dist(s.p, s.q)
    if len == 0.0
        return 0.0
    end
    return fp_dist(s.p, qr) / len
end

"""
    fp_dist_point_segment( qr, s_p, s_q )

Distance from `qr` to the segment `s_p`-`s_q`, computed from the quadratic
`a t^2 + b t + c` without constructing the nearest point. Values below
`1e-20` (squared) are reported as zero.
"""
@inline function fp_dist_point_segment(
    qr::Point{D,T},
    s_p::Point{D,T},
    s_q::Point{D,T},
) where {D,T}
    dff = s_p - qr

    a = fp_dist_sq(s_p, s_q)
    b = 2.0 * fp_dot(dff, s_q - s_p)
    c = fp_dot(dff, dff)

    local val::Float64
    if (a == 0.0)
        # A zero-length segment is a single point; -b / (2a) below would be NaN.
        val = c
    else
        t = -b / (2.0 * a)
        if (t < 0)
            val = c
        elseif (t > 1)
            val = a + b + c
        else
            val = a*t^2 + b*t + c
            if (val < 1.0e-8 * c)
                # a*t^2 + b*t cancels against c: the digits are lost when the point
                # is much closer to the segment than the segment is long. Recompute
                # from the vector to the closest point instead.
                r = dff + t * (s_q - s_p)
                val = fp_dot(r, r)
            end
        end
    end
    if (val <= 1e-20)
        return 0.0
    end
    return sqrt(val)
end

"""
    fp_dist_segment_segment( a_p, a_q, b_p, b_q )

Minimum distance between the segments `a_p`-`a_q` and `b_p`-`b_q`. Solves
the 2x2 normal equations, snaps the solution to the segments, and takes the
minimum with the four endpoint-to-segment distances. Parallel (rank one)
systems fall back to the endpoint distances.
"""
function fp_dist_segment_segment(
    a_p::Point{D,T},
    a_q::Point{D,T},
    b_p::Point{D,T},
    b_q::Point{D,T},
) where {D,T}
    # v(s,t) = (a_p - b_p)  +  s * (a_q-a_p)  +  t * (b_p-b_q)
    #        = v_1  +  s * v_2  +  t * v_3
    v_1 = a_p - b_p
    v_2 = a_q - a_p
    v_3 = b_p - b_q

    # Setting the gradient of ||v(s,t)||^2 to zero:
    #   -<v_1,v_2> = s * ||v_2||^2  +  t * <v_2,v_3>
    #   -<v_1,v_3> = s * <v_2,v_3>  +  t * ||v_3||^2
    c = [(- fp_dot(v_1, v_2)), (- fp_dot(v_1, v_3))]

    m = [
        (fp_dot(v_2, v_2)) (fp_dot(v_2, v_3));
        (fp_dot(v_2, v_3)) (fp_dot(v_3, v_3))
    ]

    rk = rank(m)
    if (rk <= 1)
        # Parallel segments, or one or both of zero length: the minimum distance is
        # realized by one of the endpoints.
        return min(
            fp_dist_point_segment(b_p, a_p, a_q),
            fp_dist_point_segment(b_q, a_p, a_q),
            fp_dist_point_segment(a_p, b_p, b_q),
            fp_dist_point_segment(a_q, b_p, b_q),
        )
    end

    b = m \ c
    # Snap solution if needed to the [0,1.0] interval....
    s = max(min(b[1], 1.0), 0.0)
    t = max(min(b[2], 1.0), 0.0)

    d::Float64 = fp_dist(fp_convex_comb(a_p, a_q, s), fp_convex_comb(b_p, b_q, t))

    return min(
        d,
        fp_dist_point_segment(b_p, a_p, a_q),
        fp_dist_point_segment(b_q, a_p, a_q),
        fp_dist_point_segment(a_p, b_p, b_q),
        fp_dist_point_segment(a_q, b_p, b_q),
    )
end

"""
    fp_geom_length

Total length of a polygonal curve, or the length of a segment.
"""
function fp_geom_length(P::PntSeq{D,T}) where {D,T}
    len::Float64 = 0.0
    for i = firstindex(P):(lastindex(P)-1)
        len += fp_dist(P[i], P[i+1])
    end

    return len
end

@inline fp_geom_length(seg::Segment{D,T}) where {D,T} = fp_dist(seg.p, seg.q);

"""
    fp_prefix_lengths

The arc length of `P` up to each vertex.
"""
function fp_prefix_lengths(P::PntSeq{D,T})::Vector{Float64} where {D,T}
    v = Vector{Float64}()
    push!(v, zero(Float64))
    len::Float64 = 0.0
    for i = firstindex(P):(lastindex(P)-1)
        len += fp_dist(P[i], P[i+1])
        push!(v, Float64(len))
    end
    return v
end

"""
    fp_match_price

The price of matching the edge p_a-p_b to the edge q_a-q_b.
"""
function fp_match_price(
    p_a::Point{N,T},
    p_b::Point{N,T},
    q_a::Point{N,T},
    q_b::Point{N,T},
) where {N,T}
    l_avg = (fp_dist(p_a, q_a) + fp_dist(p_b, q_b)) / 2.0

    return (fp_dist(p_a, p_b) + fp_dist(q_a, q_b)) * l_avg
end
