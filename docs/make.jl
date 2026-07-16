using Documenter
using Mimosa

makedocs(;
    sitename="Mimosa.jl",
    authors="Mimosa contributors",
    modules=[Mimosa],
    checkdocs=:exports,
    warnonly=[:missing_docs],
    format=Documenter.HTML(;
        canonical="https://ubercomrade.github.io/Mimosa.jl/stable/",
        prettyurls=get(ENV, "CI", nothing) == "true",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Method" => "method.md",
        "Quick Start" => "quickstart.md",
        "Julia API" => "api.md",
        "CLI" => "cli.md",
        "Supported Models" => "models.md",
        "Feature Matrix" => "feature_matrix.md",
        "Data Layout" => "data_layout.md",
        "Numerical Compatibility" => "numerical_compatibility.md",
        "Reproducibility" => "reproducibility.md",
        "Storage Format" => "storage.md",
        "Security" => "security.md",
        "Historical Python Migration" => "migration.md",
        "Extending Mimosa" => "extending.md",
        "MotifHORDE Contract" => "downstream_contract.md",
        "Architecture" => "architecture.md",
        "Release" => "release.md",
    ],
)

if get(ENV, "DOCUMENTER_DEPLOY", "false") == "true"
    deploydocs(; repo="github.com/ubercomrade/Mimosa.jl.git", devbranch="main")
end
