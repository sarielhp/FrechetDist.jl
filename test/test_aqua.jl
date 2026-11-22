module TestAqua

# Package quality checks: method ambiguities, unbound type parameters, undefined
# exports, stale dependencies, compat bounds, type piracy, persistent tasks.
using Aqua
using FrechetDist
using Test

@testset "Aqua" begin
    Aqua.test_all(FrechetDist)
end

end # module
