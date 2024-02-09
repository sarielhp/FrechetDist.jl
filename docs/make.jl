push!(LOAD_PATH,"../src/")
using FrechetDist
using Documenter

makedocs(
         sitename = "FrechetDist.jl",
         modules  = [FrechetDist],
         pages=[
                "Home" => "index.md"
         ],
    warnonly=true
)

deploydocs(;
    repo="github.com/sarielhp/FrechetDist.jl",
               )
 
