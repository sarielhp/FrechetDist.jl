#! julia

push!(LOAD_PATH, pwd()*"/src/")

using FrechetDist, FrechetDist.cg.polygon, FrechetDist.cg.point, Plots


P = Polygon2F() |> (2.0,1.0) |> (2.0, 3.4);

println( P );
