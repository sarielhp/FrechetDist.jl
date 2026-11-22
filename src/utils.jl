# Originally contributed by S. Har-Peled
# under MIT License

### Functions that deal with floating-point issues.


function fp_equal(a, b)::Bool
    return a == b
end

# Equal up to rounding. Used on arc-length positions along a curve (parameterization
# breakpoints), where two computations of the same point differ by a few ulps; the
# tolerance must stay far below the size of real features, which are small relative to
# the position on a long curve (1e-6 merged features 2 units apart at position 1e6).
function fp_equal(a::Float64, b::Float64)::Bool
    if a == b
        return true
    end
    return abs(a - b) <= (1.0e-12 * (abs(a) + abs(b)))
end
