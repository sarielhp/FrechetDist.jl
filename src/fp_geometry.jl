# Geometric primitives whose floating-point behavior the Frechet code relies
# on. Lengths, prefix lengths and norms come from BasicCompGeometry and
# LinearAlgebra, which give the same results. The functions here stay local
# because the BasicCompGeometry versions (dist, dist_sq, dot, convex_comb,
# nn_point, convex_coef, dist_point_segment, dist_segment_segment) do not give
# the same results inside the Frechet code: swapping any of them changes the
# snapshots (test/snapshots/compare.rb), and the last-bit differences change
# refinement counts and tie-breaking in the priority queues. This holds even
# for `fp_dot`, `fp_dist` and `fp_dist_sq`, which agree with BCG on isolated
# random inputs. Beyond rounding, `fp_convex_comb` returns the endpoints
# exactly at t = 0 and t = 1, `fp_dist_point_segment` guards against
# cancellation and reports tiny values as zero, and `fp_dist_segment_segment`
# uses a rank test.
#
# The kernels take static points (`Point{D,T}`, an `SVector`), whose loops the compiler
# unrolls for each D. That is fastest up to D = 16. From D = 128 on, plain loops over
# `Vector` points run 2-4.5 times faster and compile once for all D, while the static
# code takes seconds to compile at D = 1000 and minutes at D = 6000. Creating a point
# (`fp_convex_comb`) is the exception: as a `Vector` it allocates, which keeps static
# points ahead up to D = 512. Measured by example_higher/point_types.jl; see
# example_higher/README.md. The loop versions for `Vector` points are in
# fp_geometry_vec.jl.

@inline function fp_dot(p::Point{D,T}, q::Point{D,T}) where {D,T}
    s = zero(T)
    # In bounds: a Point{D,T} has exactly D coordinates.
    @inbounds for i in 1:D
        s += p[i] * q[i]
    end
    return s
end

@inline function fp_dist_sq(p::Point{D,T}, q::Point{D,T})::T where {D,T}
    sum = zero(T)
    # In bounds: a Point{D,T} has exactly D coordinates.
    @inbounds for i in 1:D
        @fastmath sum += (p[i] - q[i])^2 # julia-audit: allow-fastmath (bit-identical to the original numerics)
    end
    return sum
end

@inline function fp_dist(p::Point{D,T}, q::Point{D,T}) where {D,T}
    sum = 0.0
    # In bounds: a Point{D,T} has exactly D coordinates.
    @inbounds for i in 1:D
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
    (t == 0.0) && return p  # exact: the endpoints are returned unrounded
    (t == 1.0) && return q  # exact, as above

    o = MVector{D,T}(undef)
    # In bounds: a Point{D,T} has exactly D coordinates.
    @inbounds for i in 1:D
        @fastmath o[i] = p[i] + (q[i] - p[i]) * t # julia-audit: allow-fastmath (bit-identical to the original numerics)
        @assert(! isnan(o[i]))
    end

    return Point{D,T}(o)
end

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
    if (a == 0.0)  # exact: a zero-length segment; -b / (2a) needs a != 0
        return s_p, zero(T)
    end
    b = 2.0 * fp_dot(s_p - qr, s_vec)
    t::T = -b / (2.0 * a)

    if (t <= 0.0)
        return s_p, zero(T)
    end
    if (t >= 1.0)
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
    if (a == 0.0)  # exact: a zero-length segment; -b / (2a) needs a != 0
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

"""
    fp_convex_coef( s_p, s_q, qr )

Coefficient `d(s_p, qr) / |s_p s_q|`, or 0 for a degenerate segment.
"""
@inline function fp_convex_coef(
    s_p::AbstractVector{T},
    s_q::AbstractVector{T},
    qr::AbstractVector{T},
) where {T}
    len = fp_dist(s_p, s_q)
    if len == 0.0  # exact: a zero-length segment; avoids 0 / 0
        return 0.0
    end
    return fp_dist(s_p, qr) / len
end

"Squared relative threshold below which `fp_dist_point_segment` reports zero."
const FP_ZERO_SQ = (16 * eps(Float64))^2

"""
    fp_dist_point_segment( qr, s_p, s_q )

Distance from `qr` to the segment `s_p`-`s_q`, computed from the quadratic
`a t^2 + b t + c` without constructing the nearest point. Values below
`FP_ZERO_SQ max(a, c)` (squared rounding noise relative to the scale of the inputs)
are reported as zero.
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
            # |s_q - qr|^2 directly: a + b + c cancels when qr is near the end of a
            # long segment.
            val = fp_dist_sq(qr, s_q)
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
    # Report rounding noise as zero. The threshold is relative to the squared scale,
    # so it must be squared rounding noise ((16 eps)^2, about 1e-29), not 1e-20: that
    # zeroed every distance below 1e-10 of the segment length.
    if (val <= FP_ZERO_SQ * max(a, c))
        return 0.0
    end
    return sqrt(val)
end

"""
    fp_solve_sym2( m11, m12, m22, c1, c2 )

The solution of `[m11 m12; m12 m22] \\ [c1, c2]`, computed with the same operations
as the dense solve (a diagonal solve when `m12 == 0`, else LU with partial
pivoting that swaps only on a strictly larger pivot, and reciprocal-pivot
scaling), so the result is bit-identical, without allocating.
"""
@inline function fp_solve_sym2(
    m11::Float64,
    m12::Float64,
    m22::Float64,
    c1::Float64,
    c2::Float64,
)
    m12 == 0 && return c1 / m11, c2 / m22
    if abs(m12) > abs(m11)
        a11, a12, a21, a22, b1, b2 = m12, m22, m11, m12, c2, c1
    else
        a11, a12, a21, a22, b1, b2 = m11, m12, m12, m22, c1, c2
    end
    l = a21 * (1.0 / a11)
    u22 = a22 - l * a12
    y2 = b2 - l * b1
    x2 = y2 / u22
    x1 = (b1 - a12 * x2) / a11
    return x1, x2
end

"""
    fp_dist_segment_segment( a_p, a_q, b_p, b_q )

Minimum distance between the segments `a_p`-`a_q` and `b_p`-`b_q`. Solves
the 2x2 normal equations, snaps the solution to the segments, and takes the
minimum with the four endpoint-to-segment distances. Parallel (rank one)
systems fall back to the endpoint distances.

When `det > 1e-8 trace^2` the system certainly has rank two (`rank` would need
sigma_2 / sigma_1 <= 2 eps), and it is solved in closed form without allocating.
Otherwise (nearly parallel or degenerate segments, overflow) the original
`rank` / LAPACK path in `fp_dist_segment_segment_lapack` decides. Both give
identical results; the closed form is about 30 times faster.
"""
function fp_dist_segment_segment(
    a_p::Point{D,T},
    a_q::Point{D,T},
    b_p::Point{D,T},
    b_q::Point{D,T},
) where {D,T}
    v_1 = a_p - b_p
    v_2 = a_q - a_p
    v_3 = b_p - b_q
    m11 = fp_dot(v_2, v_2)
    m12 = fp_dot(v_2, v_3)
    m22 = fp_dot(v_3, v_3)
    tr = m11 + m22
    if !(m11 * m22 - m12 * m12 > 1.0e-8 * tr * tr)
        return fp_dist_segment_segment_lapack(a_p, a_q, b_p, b_q)
    end
    x1, x2 = fp_solve_sym2(m11, m12, m22, -fp_dot(v_1, v_2), -fp_dot(v_1, v_3))
    s = max(min(x1, 1.0), 0.0)
    t = max(min(x2, 1.0), 0.0)
    d::Float64 = fp_dist(fp_convex_comb(a_p, a_q, s), fp_convex_comb(b_p, b_q, t))
    return min(
        d,
        fp_dist_point_segment(b_p, a_p, a_q),
        fp_dist_point_segment(b_q, a_p, a_q),
        fp_dist_point_segment(a_p, b_p, b_q),
        fp_dist_point_segment(a_q, b_p, b_q),
    )
end

"The original `fp_dist_segment_segment`: dense 2x2 system, `rank`, and a dense solve."
function fp_dist_segment_segment_lapack(
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
    fp_match_price

The price of matching the edge p_a-p_b to the edge q_a-q_b.
"""
function fp_match_price(
    p_a::AbstractVector{T},
    p_b::AbstractVector{T},
    q_a::AbstractVector{T},
    q_b::AbstractVector{T},
) where {T}
    l_avg = (fp_dist(p_a, q_a) + fp_dist(p_b, q_b)) / 2.0

    return (fp_dist(p_a, p_b) + fp_dist(q_a, q_b)) * l_avg
end
