# Usage: julia compare.jl A.jls B.jls   (exit 0 when every entry is identical)
using Serialization
a = deserialize(ARGS[1]);
b = deserialize(ARGS[2])
same(x, y) = x == y || (x isa Number && y isa Number && isnan(x) && isnan(y))
function eq(x, y)
    x isa Tuple &&
        y isa Tuple &&
        return length(x) == length(y) && all(eq(i, j) for (i, j) in zip(x, y))
    x isa AbstractArray &&
        y isa AbstractArray &&
        return size(x) == size(y) && all(eq(i, j) for (i, j) in zip(x, y))
    return same(x, y)
end
bad = [k for k in sort(collect(keys(a))) if !haskey(b, k) || !eq(a[k], b[k])]
println(length(a), " entries; mismatches: ", length(bad))
foreach(k -> println("  ", k), first(bad, 25))
exit(isempty(bad) ? 0 : 1)
