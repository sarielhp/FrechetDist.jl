#! julia

push!(LOAD_PATH, pwd()*"/src/")

using FrechetDist
using BasicCompGeometry

p = point(Float64(2), Float64(4));
println( p );

pa = point(Float64(2), Float64(4), Float64( 9 ) )
println( pa );

