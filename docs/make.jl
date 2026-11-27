push!(LOAD_PATH, "../src/")
push!(LOAD_PATH, "src/")
using FrechetDist
using Documenter


makedocs(
    sitename = "FrechetDist.jl",
    modules = [FrechetDist],
    pages = [
        "Home" => "index.md",
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
