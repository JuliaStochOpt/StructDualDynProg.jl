using Documenter, StochasticDualDynamicProgramming

makedocs()

deploydocs(
    deps = Deps.pip("mkdocs", "python-markdown-math"),
    repo = "github.com/blegat/StochasticDualDynamicProgramming.jl.git",
    julia  = "release"
)
