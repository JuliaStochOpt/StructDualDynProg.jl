using Documenter, StructDualDynProg

makedocs(
    sitename = "StructDualDynProg",
    # See https://github.com/JuliaDocs/Documenter.jl/issues/868
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages = [
        "Home" => "index.md",
        "Quick Start" => "quickstart.md",
        "Tutorial" => "tutorial.md",
    ],
    # The following ensures that we only include the docstrings from
    # this module for functions define in Base that we overwrite.
    modules = [StructDualDynProg],
)

deploydocs(
    repo   = "github.com/JuliaStochOpt/StructDualDynProg.jl.git",
    push_preview = true,
)
