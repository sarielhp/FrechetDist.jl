# The sweep distance between two segments: the integral of the leash length while
# a point moves along one segment and its partner along the other, both at
# proportional speeds.
#
# With the leash v(s) = vs + s w, s in [0, 1] (vs = q_a - p_a, w = (q_b - p_b) - vs),
#
#     |v(s)|^2 = a (s - t0)^2 + h^2,   a = |w|^2,  t0 = -<vs, w> / a,  h = |v(t0)|,
#
# so the integral is sqrt(a) * the integral of sqrt(u^2 + k^2), k = h / sqrt(a), over
# u in [-t0, 1 - t0], which has a closed form. Computing h from the vector v(t0)
# (rather than from the discriminant 4ac - b^2) and evaluating the closed form
# without subtracting nearly equal terms keeps the result accurate to a few ulps at
# every scale.

"""
    sqrt_quad_integral( u1, u2, k )

The integral of `sqrt(u^2 + k^2)` over `u` in `[u1, u2]`, for `k >= 0` and an
interval of width one (`u2 - u1 == 1` is used exactly). The antiderivative is
`(u s + k^2 asinh(u / k)) / 2` with `s = sqrt(u^2 + k^2)`; it is odd, so an interval
around 0 adds two positive terms, and an interval on one side is rewritten so that
no two nearly equal terms are subtracted.
"""
function sqrt_quad_integral(u1::Float64, u2::Float64, k::Float64)
    k2 = k * k
    if u1 < 0.0 < u2
        F(u) = (u * sqrt(u * u + k2) + (k2 == 0.0 ? 0.0 : k2 * asinh(u / k))) / 2
        return F(u2) - F(u1)
    end
    if u2 <= 0.0                 # mirror to 0 <= u1 < u2
        u1, u2 = -u2, -u1
    end
    s1 = sqrt(u1 * u1 + k2)
    s2 = sqrt(u2 * u2 + k2)
    # u2 s2 - u1 s1, using u2 - u1 = 1.
    d1 = (u2 + u1) * (u2 * u2 + u1 * u1 + k2) / (u2 * s2 + u1 * s1)
    # k^2 (asinh(u2 / k) - asinh(u1 / k)) = k^2 log((u2 + s2) / (u1 + s1)).
    d2 = (k2 == 0.0) ? 0.0 : k2 * log1p((1.0 + (u2 + u1) / (s2 + s1)) / (u1 + s1))
    return (d1 + d2) / 2
end

"""
    SweepDist_segs_p( p_a, p_b, q_a, q_b )

The oriented sweep distance (continuous dynamic time warping) between the oriented
segments `p_a p_b` and `q_a q_b`: the integral, over the arc length of `p_a p_b`, of
the distance between the point on `p_a p_b` and the corresponding point (at the same
fraction of the way) on `q_a q_b`.
"""
function SweepDist_segs_p(
    p_a::Point{D,T},
    p_b::Point{D,T},
    q_a::Point{D,T},
    q_b::Point{D,T},
) where {D,T}
    L = fp_dist(p_a, p_b)
    L == 0.0 && return 0.0

    vs = q_a - p_a
    w = (q_b - p_b) - vs
    a = fp_dot(w, w)
    a == 0.0 && return L * norm(vs)    # a translation (or identical segments)

    t0 = -fp_dot(vs, w) / a            # parameter of the shortest leash
    h = norm(vs + t0 * w)              # its length
    ra = sqrt(a)
    return L * ra * sqrt_quad_integral(-t0, 1.0 - t0, h / ra)
end

"""
    SweepDist_segs( p_a, p_b, q_a, q_b )

The sweep distance (continuous dynamic time warping) between the segments
``p_a p_b`` and ``q_a q_b``: the sum of the two oriented integrals of the distance.
"""
function SweepDist_segs(
    p_a::Point{D,T},
    p_b::Point{D,T},
    q_a::Point{D,T},
    q_b::Point{D,T},
) where {D,T}
    lp = SweepDist_segs_p(p_a, p_b, q_a, q_b)
    lq = SweepDist_segs_p(q_a, q_b, p_a, p_b)

    return lp + lq
end
