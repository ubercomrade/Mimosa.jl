using PackageCompiler

const ROOT = normpath(joinpath(@__DIR__, ".."))
const OUTPUT = abspath(get(ENV, "MIMOSA_APP_OUTPUT", joinpath(ROOT, "dist", "Mimosa")))

mkpath(dirname(OUTPUT))

create_app(
    ROOT,
    OUTPUT;
    executables=["mimosa" => "julia_main"],
    precompile_execution_file=joinpath(@__DIR__, "precompile_app.jl"),
    force=true,
)
