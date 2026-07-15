function executable_path(bundle::AbstractString)
    name = Sys.iswindows() ? "mimosa.exe" : "mimosa"
    return joinpath(abspath(bundle), "bin", name)
end

length(ARGS) == 1 || error("usage: smoke_app.jl <compiled-app-directory>")

executable = executable_path(only(ARGS))
isfile(executable) || error("compiled executable not found: $executable")

run(`$executable --version`)
run(`$executable --help`)
