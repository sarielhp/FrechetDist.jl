# The kernels of fp_geometry.jl for points whose length is known only at run time:
# `Vector` points, as in `GPntSeq{AnyDim,Float64,Vector{Float64}}`. The static methods
# (`Point{D,T}`) are more specific, so static curves never reach these.
#
# Each kernel checks once that its points have the same length and then loops in
# bounds, accumulating the dot products it needs in one pass, without temporaries.
# Only `fp_convex_comb` allocates (it creates the new point). The results agree with
# the static kernels up to rounding, not bit for bit: the sums are reassociated
# (`@simd`). Above about D = 128 these run faster than the static kernels and compile
# once for every dimension (example_higher/point_types.jl).

@noinline throw_point_lengths(m, n) =
    throw(DimensionMismatch("points with $m and $n coordinates"))

@inline function vec_length(p::AbstractVector, q::AbstractVector)
    Base.require_one_based_indexing(p, q)
    n = length(p)
    length(q) == n || throw_point_lengths(n, length(q))
    return n
end

@inline function vec_length(p::AbstractVector, q::AbstractVector, r::AbstractVector)
    n = vec_length(p, q)
    vec_length(p, r)
    return n
end

@inline function fp_dot(p::AbstractVector{T}, q::AbstractVector{T}) where {T}
    n = vec_length(p, q)
    s = zero(T)
    # In bounds: p and q have n coordinates, indexed from 1 (vec_length).
    @inbounds @simd for i in 1:n
        s += p[i] * q[i]
    end
    return s
end

@inline function fp_dist_sq(p::AbstractVector{T}, q::AbstractVector{T})::T where {T}
    n = vec_length(p, q)
    s = zero(T)
    # In bounds: as in fp_dot.
    @inbounds @simd for i in 1:n
        s += (p[i] - q[i])^2
    end
    return s
end

@inline fp_dist(p::AbstractVector{T}, q::AbstractVector{T}) where {T} =
    sqrt(fp_dist_sq(p, q))

@inline function fp_convex_comb(
    p::AbstractVector{T},
    q::AbstractVector{T},
    t::Float64,
) where {T}
    (t == 0.0) && return p  # exact: the endpoints are returned unrounded (and shared)
    (t == 1.0) && return q  # exact, as above
    n = vec_length(p, q)
    o = Vector{T}(undef, n)
    # In bounds: p, q and o have n coordinates.
    @inbounds @simd for i in 1:n
        o[i] = p[i] + (q[i] - p[i]) * t
    end
    @assert(!any(isnan, o))
    return o
end

# Coordinate i of fp_convex_comb(p, q, t), without creating the point.
@inline function vec_comb_at(p, q, t::Float64, i::Int)
    (t == 0.0) && return p[i]
    (t == 1.0) && return q[i]
    return p[i] + (q[i] - p[i]) * t
end

# <p - r, v>, without the temporary p - r.
@inline function vec_dot_diff(
    p::AbstractVector{T},
    r::AbstractVector{T},
    v::AbstractVector{T},
) where {T}
    n = vec_length(p, r, v)
    s = zero(T)
    # In bounds: as in fp_dot.
    @inbounds @simd for i in 1:n
        s += (p[i] - r[i]) * v[i]
    end
    return s
end

# <s_p - qr, s_q - s_p>, without temporaries.
@inline function vec_dot_seg(
    s_p::AbstractVector{T},
    s_q::AbstractVector{T},
    qr::AbstractVector{T},
) where {T}
    n = vec_length(s_p, s_q, qr)
    s = zero(T)
    # In bounds: as in fp_dot.
    @inbounds @simd for i in 1:n
        s += (s_p[i] - qr[i]) * (s_q[i] - s_p[i])
    end
    return s
end

@inline function fp_nn_point_ext_ext(
    s_p::AbstractVector{T},
    s_q::AbstractVector{T},
    qr::AbstractVector{T},
    seg_len_sq::T,
    s_vec::AbstractVector{T},
) where {T}
    a = seg_len_sq
    (a == 0.0) && return s_p, zero(T)  # exact: a zero-length segment
    t = T(-2.0 * vec_dot_diff(s_p, qr, s_vec) / (2.0 * a))
    (t <= 0.0) && return s_p, zero(T)
    (t >= 1.0) && return s_q, one(T)
    return fp_convex_comb(s_p, s_q, t), t
end

@inline function fp_nn_point_ext(
    s_p::AbstractVector{T},
    s_q::AbstractVector{T},
    qr::AbstractVector{T},
) where {T}
    a = fp_dist_sq(s_p, s_q)
    (a == 0.0) && return s_p, zero(T)  # exact: a zero-length segment
    t = T(-2.0 * vec_dot_seg(s_p, s_q, qr) / (2.0 * a))
    (t <= 0.0) && return s_p, zero(T)
    (t >= 1.0) && return s_q, one(T)
    return fp_convex_comb(s_p, s_q, t), t
end

@inline function fp_nn_point(
    s_p::AbstractVector{T},
    s_q::AbstractVector{T},
    qr::AbstractVector{T},
) where {T}
    return first(fp_nn_point_ext(s_p, s_q, qr))
end

@inline function fp_dist_point_segment(
    qr::AbstractVector{T},
    s_p::AbstractVector{T},
    s_q::AbstractVector{T},
) where {T}
    n = vec_length(qr, s_p, s_q)
    a = b2 = c = zero(T)
    # In bounds: as in fp_dot. dff = s_p - qr, v = s_q - s_p.
    @inbounds @simd for i in 1:n
        dff = s_p[i] - qr[i]
        v = s_q[i] - s_p[i]
        a += v * v
        b2 += dff * v
        c += dff * dff
    end
    b = 2.0 * b2
    local val::Float64
    if (a == 0.0)
        val = c
    else
        t = -b / (2.0 * a)
        if (t < 0)
            val = c
        elseif (t > 1)
            val = fp_dist_sq(qr, s_q)
        else
            val = a*t^2 + b*t + c
            if (val < 1.0e-8 * c)
                val = zero(T)
                # In bounds: as above.
                @inbounds @simd for i in 1:n
                    val += (s_p[i] - qr[i] + t * (s_q[i] - s_p[i]))^2
                end
            end
        end
    end
    # Rounding noise is reported as zero, as in the static kernel.
    (val <= FP_ZERO_SQ * max(a, c)) && return 0.0
    return sqrt(val)
end

function fp_dist_segment_segment(
    a_p::AbstractVector{T},
    a_q::AbstractVector{T},
    b_p::AbstractVector{T},
    b_q::AbstractVector{T},
) where {T}
    n = vec_length(a_p, a_q, b_p)
    vec_length(a_p, b_q)
    m11 = m12 = m22 = c1 = c2 = zero(T)
    # In bounds: as in fp_dot. v_1 = a_p - b_p, v_2 = a_q - a_p, v_3 = b_p - b_q.
    @inbounds @simd for i in 1:n
        v_1 = a_p[i] - b_p[i]
        v_2 = a_q[i] - a_p[i]
        v_3 = b_p[i] - b_q[i]
        m11 += v_2 * v_2
        m12 += v_2 * v_3
        m22 += v_3 * v_3
        c1 += v_1 * v_2
        c2 += v_1 * v_3
    end
    ends = min(
        fp_dist_point_segment(b_p, a_p, a_q),
        fp_dist_point_segment(b_q, a_p, a_q),
        fp_dist_point_segment(a_p, b_p, b_q),
        fp_dist_point_segment(a_q, b_p, b_q),
    )
    tr = m11 + m22
    if m11 * m22 - m12 * m12 > 1.0e-8 * tr * tr
        x1, x2 = fp_solve_sym2(m11, m12, m22, -c1, -c2)
    else
        # Nearly parallel or degenerate: the rank test of the static kernel.
        m = [m11 m12; m12 m22]
        rank(m) <= 1 && return ends
        x1, x2 = m \ [-c1, -c2]
    end
    s = max(min(x1, 1.0), 0.0)
    t = max(min(x2, 1.0), 0.0)
    d = zero(T)
    # In bounds: as above.
    @inbounds for i in 1:n
        d += (vec_comb_at(a_p, a_q, s, i) - vec_comb_at(b_p, b_q, t, i))^2
    end
    return min(sqrt(d), ends)
end
