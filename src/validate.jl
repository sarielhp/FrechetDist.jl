# Input validation for the public entry points.
#
# Every exported function that takes curves checks them up front and throws an
# ArgumentError that names the function, instead of failing later (or silently
# returning a value) with an AssertionError, BoundsError, KeyError, or NaN.

"Largest coordinate magnitude accepted: squared distances must not overflow."
const MAX_COORDINATE = 1.0e150

"""
Largest number of vertices: event ids pack a vertex index into 30 bits, and the search
forms ids for index `n + 1` (one past the last vertex), which must still fit.
"""
const MAX_VERTICES = 2^30 - 2

"""
Largest number of refinement rounds the iterative algorithms (`frechet_c_compute`,
`frechet_c_approx`, `frechet_mono_via_refinement`) may run before giving up with an
error instead of looping forever. At most 17 rounds were observed in about 10,000 runs.
It is a `Ref` so that tests can lower it.
"""
const MAX_REFINEMENT_ROUNDS = Ref(1000)

function check_rounds(rounds::Int, who::AbstractString)
    rounds > MAX_REFINEMENT_ROUNDS[] &&
        error("$who: no convergence after $(MAX_REFINEMENT_ROUNDS[]) refinement rounds")
    return nothing
end

"""
    check_curve( P, who, name; min_vertices = 2 )

Throws an `ArgumentError` mentioning `who` unless the curve `P` has `Float64`
coordinates, between `min_vertices` and `MAX_VERTICES` vertices, and all
coordinates finite with magnitude at most `MAX_COORDINATE`.
"""
function check_curve(
    P::GPntSeq{D,T,PT},
    who::AbstractString,
    name::AbstractString;
    min_vertices::Int = 2,
) where {D,T,PT<:AbstractVector{T}}
    check_coordinate_type(P, who, name)
    check_point_lengths(P, who, name)
    n = cardin(P)
    n >= min_vertices || throw(
        ArgumentError("$who: $name has $n vertices; at least $min_vertices are required"),
    )
    n <= MAX_VERTICES || throw(ArgumentError("$who: $name has too many vertices ($n)"))
    for p in P.pnts, x in p
        isfinite(x) || throw(ArgumentError("$who: $name has a non-finite coordinate ($x)"))
        abs(x) <= MAX_COORDINATE || throw(
            ArgumentError(
                "$who: a coordinate of $name ($x) exceeds $MAX_COORDINATE in magnitude",
            ),
        )
    end
    return nothing
end

"""
    check_coordinate_type( P, who, name = "the curve" )

Throws an `ArgumentError` mentioning `who` unless the coordinates of `P` are `Float64`.
Constant time: for entry points called once per edge, where `check_curve` would be
quadratic overall.
"""
function check_coordinate_type(
    ::GPntSeq{D,T,PT},
    who::AbstractString,
    name::AbstractString = "the curve",
) where {D,T,PT<:AbstractVector{T}}
    T === Float64 ||
        throw(ArgumentError("$who: the coordinates of $name must be Float64 (got $T)"))
    # The kernels create static points from static ones and `Vector`s from the rest,
    # and push them into the curve, so only these two point types work.
    (PT === Point{D,T} || PT === Vector{T}) || throw(
        ArgumentError(
            "$who: $name has points of type $PT; use Point{D,$T} (static) or " *
            "Vector{$T} points, or a MatVecPntSeq, which is converted",
        ),
    )
    # Events and morphings take their dimension parameter from the point type.
    point_dim(PT) === D || throw(
        ArgumentError(
            "$who: $name has dimension parameter $D but points of type $PT; for " *
            "points of run-time length use GPntSeq{AnyDim,$T,$PT}(P.pnts), which " *
            "shares the points",
        ),
    )
    return nothing
end

"Throws unless all points of `P` have the same number of coordinates."
check_point_lengths(
    ::GPntSeq{D,T,PT},
    ::AbstractString,
    ::AbstractString,
) where {D,T,PT<:StaticVector} = nothing
function check_point_lengths(P::GPntSeq, who::AbstractString, name::AbstractString)
    isempty(P.pnts) && return nothing
    d = length(first(P.pnts))
    for p in P.pnts
        length(p) == d || throw(
            ArgumentError("$who: $name has points with $d and $(length(p)) coordinates"),
        )
    end
    return nothing
end

function check_curves(P::GPntSeq, Q::GPntSeq, who::AbstractString; min_vertices::Int = 2)
    check_curve(P, who, "the first curve"; min_vertices)
    check_curve(Q, who, "the second curve"; min_vertices)
    # Equal by type for static points; checked here for points of run-time length.
    (isempty(P) || isempty(Q) || dim(P) == dim(Q)) || throw(
        ArgumentError(
            "$who: the two curves must have the same dimension (got $(dim(P)) and $(dim(Q)))",
        ),
    )
    return nothing
end

"""
    check_approx( approx, who; exceed_one = true )

Returns `approx` as a `Float64`; throws unless it is finite, and larger than 1
when `exceed_one`.
"""
function check_approx(approx::Real, who::AbstractString; exceed_one::Bool = true)
    a = Float64(approx)
    isfinite(a) || throw(ArgumentError("$who: approx must be finite (got $approx)"))
    (exceed_one && !(a > 1.0)) &&
        throw(ArgumentError("$who: approx must be larger than 1 (got $approx)"))
    return a
end

"Returns `tolerance` as a `Float64`; throws unless it is finite and non-negative."
function check_tolerance(tolerance::Real, who::AbstractString)
    t = Float64(tolerance)
    (isfinite(t) && t >= 0.0) || throw(
        ArgumentError("$who: tolerance must be finite and non-negative (got $tolerance)"),
    )
    return t
end

"""
    check_offsets( offs, who, what )

Throws unless every entry of `offs` is finite and non-negative.
"""
function check_offsets(offs::Vector{Float64}, who::AbstractString, what::AbstractString)
    for x in offs
        (isfinite(x) && x >= 0) || throw(
            ArgumentError(
                "$who: the entries of $what must be finite and non-negative (got $x)",
            ),
        )
    end
    return nothing
end

"Throws unless `x` (a width, bound, or radius) is finite and non-negative."
function check_nonnegative(x::Real, who::AbstractString, what::AbstractString)
    (isfinite(x) && x >= 0) ||
        throw(ArgumentError("$who: $what must be finite and non-negative (got $x)"))
    return nothing
end

"""
    check_quadratic_size( n, m, bytes_per_cell, who )

Throws an `ArgumentError` when an algorithm that allocates `bytes_per_cell` bytes for each
of the `n * m` vertex pairs needs more than the total memory of the machine, instead of
failing with an out-of-memory error after a long allocation.
"""
function check_quadratic_size(n::Int, m::Int, bytes_per_cell::Real, who::AbstractString)
    needed = Float64(n) * Float64(m) * bytes_per_cell
    available = Float64(Sys.total_memory())
    needed <= available || throw(
        ArgumentError(
            "$who needs about $(round(needed / 2^30; digits = 1)) GiB for curves of " *
            "$n and $m vertices, but this machine has $(round(available / 2^30; digits = 1)) GiB",
        ),
    )
    return nothing
end

"""
    check_indices( indices, n, who )

Throws an `ArgumentError` unless `indices` are vertex indices in `1:n` in non-decreasing
order.
"""
function check_indices(indices::Vector{Int64}, n::Int, who::AbstractString)
    for (k, i) in enumerate(indices)
        (1 <= i <= n) || throw(ArgumentError("$who: vertex index $i is outside 1:$n"))
        (k > 1 && i < indices[k-1]) &&
            throw(ArgumentError("$who: vertex indices must be in non-decreasing order"))
    end
    return nothing
end

"Whether the points of `P` are static vectors (their length is part of the type)."
static_points(P::AbsPntSeq) = eltype(P.pnts) <: StaticVector

"`P` with static points, as a `PntSeq` (a copy unless it is one already)."
static_curve(P::GPntSeq) = P
static_curve(P::AbsPntSeq{D,T}) where {D,T} = PntSeq{D,T}(collect(P.pnts))

"`P` with `Vector` points of run-time length (a copy of the coordinates unless it is one already)."
vector_curve(P::GPntSeq{AnyDim,T,Vector{T}}) where {T} = P
vector_curve(P::AbsPntSeq{D,T}) where {D,T} =
    GPntSeq{AnyDim,T,Vector{T}}([Vector{T}(p) for p in P.pnts])

"""
    common_curves( P, Q )

`P` and `Q` as two curves of one type that the algorithms take, or `nothing` when the
coordinate types differ or the points are static of different dimensions. Static
points stay static. Otherwise (a `MatVecPntSeq`, a matrix column view, or static and
`Vector` points mixed) both are copied into `GPntSeq{AnyDim,T,Vector{T}}`.
"""
function common_curves(P::AbsPntSeq{D1,T}, Q::AbsPntSeq{D2,T}) where {D1,D2,T}
    if static_points(P) && static_points(Q)
        D1 == D2 || return nothing
        return static_curve(P), static_curve(Q)
    end
    return vector_curve(P), vector_curve(Q)
end
common_curves(::AbsPntSeq, ::AbsPntSeq) = nothing

"`P` as a curve of a type the algorithms take (see `common_curves`)."
single_curve(P::AbsPntSeq) = static_points(P) ? static_curve(P) : vector_curve(P)

# Curves of other representations (a `MatVecPntSeq`, mixed point types) are converted
# once, at the entry, by `common_curves`. Curves of different dimension or coordinate
# type do not match the typed methods; the error says so instead of reporting a
# MethodError listing every method. Curves of a type the algorithms take, with wrong
# further arguments, give a plain MethodError.
for f in (
    :frechet_c_compute,
    :frechet_c_decide_exact,
    :frechet_c_compute_exact,
    :frechet_c_approx,
    :frechet_d_compute,
    :frechet_d_r_compute,
    :DTW_d_compute,
    :frechet_ve_r_compute,
    :FEVER_compute,
    :SweepDist_compute,
    :SweepDist_lb_compute,
    :frechet_mono_via_refinement,
    :frechet_dist_upper_bound,
)
    @eval function $f(P::AbsPntSeq, Q::AbsPntSeq, args...; kwargs...)
        C = common_curves(P, Q)
        if C !== nothing
            # Converted curves are of one type the algorithms take, so a second pass
            # here (wrong further arguments) ends with the MethodError.
            typeof(C) === typeof((P, Q)) && throw(MethodError($f, (P, Q, args...)))
            return $f(C[1], C[2], args...; kwargs...)
        end
        throw(
            ArgumentError(
                $(string(f)) *
                ": the two curves must have the same dimension and coordinate type " *
                "(got $(typeof(P)) and $(typeof(Q)))",
            ),
        )
    end
end

# Single-curve entries: other representations are converted as above.
for f in (
    :frechet_simplify_to_width,
    :frechet_simplify_to_cardin,
    :frechet_simplify_w_exp,
    :frechet_width_approx,
)
    @eval function $f(P::AbsPntSeq, args...; kwargs...)
        C = single_curve(P)
        C === P && throw(MethodError($f, (P, args...)))
        return $f(C, args...; kwargs...)
    end
end
