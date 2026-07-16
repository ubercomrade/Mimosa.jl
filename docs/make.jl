using Documenter
using Mimosa

makedocs(;
    sitename="Mimosa.jl",
    authors="Mimosa contributors",
    modules=[Mimosa],
    # The guide documents task-oriented entry points rather than every export.
    checkdocs=:none,
    format=Documenter.HTML(;
        canonical="https://ubercomrade.github.io/Mimosa.jl/stable/",
        prettyurls=get(ENV, "CI", nothing) == "true",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Quick Start" => "quickstart.md",
        "Julia API Guide" => "api.md",
        "CLI" => "cli.md",
        "Supported Models" => "models.md",
        "Method" => "method.md",
        "Data Layout" => "data_layout.md",
        "Storage Format" => "storage.md",
    ],
)

if get(ENV, "DOCUMENTER_DEPLOY", "false") == "true"
    deploydocs(; repo="github.com/ubercomrade/Mimosa.jl.git", devbranch="main")
end
