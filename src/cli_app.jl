# Entry point used by `Pkg.Apps` (`pkg> app add Mimosa`).
module CLIApp

using ..Mimosa

function (@main)(args::Vector{String})
    return Mimosa.main(args)
end

end # module CLIApp

# Entry point used by `PackageCompiler.create_app`.
function julia_main()::Cint
    return Cint(main(ARGS))
end
