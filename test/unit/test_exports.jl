# Test that every exported symbol has a docstring.
# Fail-closed: missing docstrings on exported names cause test failure.

using Test
using Mimosa

@testset "Exported symbols have docstrings" begin
    exported = names(Mimosa; all=false)
    @test !isempty(exported)

    # Use Docs.meta to check which symbols have docstrings
    meta = Docs.meta(Mimosa)
    documented_vars = Set{Symbol}()
    for k in keys(meta)
        push!(documented_vars, k.var)
    end

    missing_docs = String[]
    for sym in exported
        name_str = string(sym)
        startswith(name_str, "_") && continue

        if !(sym in documented_vars)
            push!(missing_docs, name_str)
        end
    end

    if !isempty(missing_docs)
        @error "Exported symbols missing docstrings" missing=missing_docs
    end
    @test isempty(missing_docs)
end
