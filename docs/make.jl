using ExcitonsBSE
using Documenter

DocMeta.setdocmeta!(ExcitonsBSE, :DocTestSetup, :(using ExcitonsBSE); recursive=true)

makedocs(;
    modules=[ExcitonsBSE],
    authors="Bruno Amorim <amorim.bac@gmail.com> and contributors",
    sitename="ExcitonsBSE.jl",
    format=Documenter.HTML(;
        canonical="https://BacAmorim.github.io/ExcitonsBSE.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/BacAmorim/ExcitonsBSE.jl",
    devbranch="main",
)
