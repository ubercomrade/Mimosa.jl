using Test
using Mimosa

# Aqua quality checks — fail-closed (Aqua is a required test dependency).
using Aqua
Aqua.test_all(Mimosa; stale_deps=false, project_extras=false)

@testset "Mimosa.jl Stage 1-7" begin
    # Unit tests
    include("unit/test_models.jl")
    include("unit/test_readers.jl")
    include("unit/test_serialization.jl")
    include("unit/test_sequences.jl")
    include("unit/test_profiles.jl")
    include("unit/test_sites.jl")
    include("unit/test_bamm.jl")
    include("unit/test_sitega.jl")
    include("unit/test_dimont.jl")
    include("unit/test_slim.jl")
    include("unit/test_gev.jl")
    include("unit/test_pvalues.jl")
    include("unit/test_relations.jl")
    include("unit/test_null_distribution.jl")
    include("unit/test_null_storage.jl")
    include("unit/test_parallel.jl")
    include("unit/test_cache.jl")
    include("unit/test_model_storage.jl")
    include("unit/test_validation.jl")
    include("unit/test_model_geometry.jl")
    include("unit/test_extending.jl")
    include("unit/test_exports.jl")
    include("unit/test_type_stability.jl")

    # Property tests
    include("properties/test_properties.jl")

    # Integration tests (CLI path)
    include("integration/test_cli.jl")
    include("integration/test_cli_subprocess.jl")

    # JET static analysis (fail-closed for type instability in hot paths)
    include("jet/test_jet.jl")
end

# Downstream contract test — runs in a separate consumer environment:
#   julia --project=test/downstream test/downstream/runtests.jl
