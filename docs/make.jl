push!(LOAD_PATH,"../src/")
push!(LOAD_PATH,"src/")
import Pkg;
Pkg.add("StaticArrays")
Pkg.add("DataStructures")

using FrechetDist
using Documenter


makedocs(
         sitename = "FrechetDist.jl",
         modules  = [FrechetDist],
         pages=[
             "Home" => "index.md",
             "cg.md",
             "frechet_discrete.md",
             "morphing.md"
         ],
    warnonly=true
)

deploydocs(;
    repo="github.com/sarielhp/FrechetDist.jl",
               )
 
