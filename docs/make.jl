using RiskBudgeting
using Documenter

DocMeta.setdocmeta!(RiskBudgeting, :DocTestSetup, :(using RiskBudgeting); recursive=true)

makedocs(;
    modules=[RiskBudgeting],
    authors="Kiar Fatah",
    repo="https://github.com/xiar-fatah/RiskBudgeting.jl/blob/{commit}{path}#{line}",
    sitename="RiskBudgeting.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://xiar-fatah.github.io/RiskBudgeting.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/xiar-fatah/RiskBudgeting.jl",
    devbranch="main",
)
