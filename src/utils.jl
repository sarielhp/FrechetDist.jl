# Originally contributed by S. Har-Peled
# under MIT License

### Functions try to deal with floating point issues.


function fp_equal(a, b)::Bool
    return a == b
end

function fp_equal(a::Float64, b::Float64)::Bool
    if a == b
        return true
    end
    return abs(a - b) <= (0.000001 * (abs(a) + abs(b)))
end

function fp_equal(p::Point{D,T}, q::Point{D,T})::Bool where {D,T}
    for i = 1:D
        if ! fp_equal(p[i], q[i])
            return false
        end
    end
    return true
end
