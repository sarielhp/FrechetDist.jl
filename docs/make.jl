push!(LOAD_PATH, "../src/")
push!(LOAD_PATH, "src/")
using FrechetDist
using Documenter


makedocs(
    sitename = "FrechetDist.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
        footer = "[FrechetDist.jl on GitHub](https://github.com/sarielhp/FrechetDist.jl). " *
                 "Powered by [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl).",
    ),
    modules = [FrechetDist],
    pages = [
        "Home" => "index.md",
        "installation.md",
        "quickstart.md",
        "algorithms.md",
        "cg.md",
        "frechet_discrete.md",
        "morphing.md",
        "simplification.md",
        "sweep.md",
        "low_level.md",
        "examples.md",
        "high_dimensional.md",
        "changelog.md",
    ],
    # An exported name missing from the manual fails the build; other problems warn.
    checkdocs = :exports,
    warnonly = Documenter.except(:missing_docs),
)

deploydocs(;
    repo = "github.com/sarielhp/FrechetDist.jl.git",
    devbranch = "main",
    forcepush = true,
    push_preview = true,
)
