#! julia

push!(LOAD_PATH, pwd()*"/src/")

using FrechetDist
using FrechetDist.cg
using FrechetDist.cg.point

p = npoint(Float64(2), Float64(4));
println( p );

pa = npoint(Float64(2), Float64(4), Float64( 9 ) )
println( pa );

