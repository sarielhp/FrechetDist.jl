# Input validation for the public entry points.
#
# Every exported function that takes curves checks them up front and throws an
# ArgumentError that names the function, instead of failing later (or silently
# returning a value) with an AssertionError, BoundsError, KeyError, or NaN.

"Largest coordinate magnitude accepted: squared distances must not overflow."
const MAX_COORDINATE = 1.0e150

"Largest number of vertices: event ids pack a vertex index into 30 bits."
const MAX_VERTICES = 2^30 - 1

"""
Largest number of refinement rounds the iterative algorithms (`frechet_c_compute`,
`frechet_c_approx`, `frechet_mono_via_refinement`) may run before giving up with an
error instead of looping forever. At most 17 rounds were observed in about 10,000 runs.
A `Ref` so that tests can lower it.
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
    P::PntSeq{D,T},
    who::AbstractString,
    name::AbstractString;
    min_vertices::Int = 2,
) where {D,T}
    T === Float64 ||
        throw(ArgumentError("$who: the coordinates of $name must be Float64 (got $T)"))
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

function check_curves(P::PntSeq, Q::PntSeq, who::AbstractString; min_vertices::Int = 2)
    check_curve(P, who, "the first curve"; min_vertices)
    check_curve(Q, who, "the second curve"; min_vertices)
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

# Curves of different dimension or coordinate type do not match the typed methods;
# say so instead of reporting a MethodError listing every method.
for f in (
    :frechet_c_compute,
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
    @eval function $f(P::PntSeq, Q::PntSeq, args...; kwargs...)
        throw(
            ArgumentError(
                $(string(f)) *
                ": the two curves must have the same dimension and coordinate type " *
                "(got $(typeof(P)) and $(typeof(Q)))",
            ),
        )
    end
end
