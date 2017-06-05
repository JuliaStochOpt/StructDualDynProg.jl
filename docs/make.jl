using Documenter, StructDualDynProg

makedocs()

deploydocs(
    deps = Deps.pip("mkdocs", "python-markdown-math"),
    repo = "github.com/blegat/StructDualDynProg.jl.git",
    julia  = "release"
)
