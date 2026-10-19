#! julia

push!(LOAD_PATH, pwd()*"/src/")

using FrechetDist, Plots
using BasicCompGeometry


P = PntSeq( [ point( 2.0, 1.0 ), point( 2.0, 3.4 ) ] );

println( P );
