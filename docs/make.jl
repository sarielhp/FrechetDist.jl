push!(LOAD_PATH,"../src/")
using FrechetDist
using Documenter

makedocs(
         sitename = "FrechetDist.jl",
         modules  = [FrechetDist],
         pages=[
                "Home" => "index.md"
               ])deploydocs(;
    repo="github.com/sarielhp/FrechetDist.jl",
)
 
