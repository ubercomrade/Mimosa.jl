using Test
using Mimosa

# Aqua quality checks — fail-closed (Aqua is a required test dependency).
using Aqua
Aqua.test_all(Mimosa; stale_deps=false, project_extras=false)

function _test_context_model_contract(Model)
    model_name = string(nameof(Model))
    model0 = Model("test0", zeros(Float32, 5, 10), 0, 10)
    @test model0.name == "test0"
    @test model0.span == 0
    @test model0.motif_length == 10
    @test size(model0) == (5, 10)
    @test length(model0) == 10
    @test eltype(model0) == Float32

    model1 = Model("test1", zeros(Float32, 25, 8), 1, 8)
    @test model1.span == 1
    @test size(model1) == (25, 8)

    model3 = Model("test3", zeros(Float32, 625, 5), 3, 5)
    @test model3.span == 3
    @test size(model3) == (625, 5)

    @test_throws MimosaError Model("bad", Matrix{Float32}(undef, 4, 3), 0, 3)
    @test_throws MimosaError Model("bad", Matrix{Float32}(undef, 24, 3), 1, 3)
    @test_throws MimosaError Model("bad", Matrix{Float32}(undef, 5, 0), 0, 0)
    @test_throws MimosaError Model("bad", Matrix{Float32}(undef, 5, 3), -1, 3)
    invalid = zeros(Float32, 5, 3)
    invalid[1, 1] = NaN32
    @test_throws MimosaError Model("bad", invalid, 0, 3)

    shown = sprint(show, Model("example", zeros(Float32, 5, 13), 0, 13))
    @test contains(shown, model_name)
    @test contains(shown, "example")
    @test contains(shown, "span=0")

    representation = ones(Float32, 5, 3)
    a = Model("x", representation, 0, 3)
    b = Model("x", copy(representation), 0, 3)
    @test a == b
    @test a != Model("y", representation, 0, 3)
    @test a != Model("x", ones(Float32, 25, 3), 1, 3)
    @test isapprox(a, b)

    bounds = Float32[
        1 2 3
        -1 -2 -3
        0.5 1.5 2.5
        0 0 0
        -1 -2 -3
    ]
    @test scorebounds(Model("bounds", bounds, 0, 3)) == (-6.0f0, 6.0f0)
end

function _test_scan_contract(model, sequence)
    n_positions = npositions(model, length(sequence))
    @test n_positions == length(sequence) - Mimosa.window_size(model) + 1

    forward = scan(model, sequence; strands=ForwardOnly())
    reverse = scan(model, sequence; strands=ReverseOnly())
    best = scan(model, sequence; strands=BestStrand())
    both = scan(model, sequence; strands=BothStrands())

    @test length(forward) == n_positions
    @test length(reverse) == n_positions
    @test all(isfinite, forward)
    @test all(isfinite, reverse)
    @test all(best .>= min.(forward, reverse))
    @test both.forward ≈ forward
    @test both.reverse ≈ reverse

    destination = similar(forward)
    scan!(destination, model, sequence; strands=ForwardOnly())
    @test destination ≈ forward
    @test isempty(scan(model, UInt8[0, 1, 2, 3]; strands=ForwardOnly()))
    @test scan(model, sequence; strands=ForwardOnly()) == forward

    sequence_copy = copy(sequence)
    scan(model, sequence; strands=ForwardOnly())
    @test sequence == sequence_copy
end

function _test_batch_scan_contract(model, sequences)
    batch = EncodedSequenceBatch([sequences[1], sequences[2], UInt8[0, 0, 0, 0]])
    expected_length = npositions(model, length(sequences[1]))

    scores = scan(model, batch; strands=ForwardOnly())
    both = scan(model, batch; strands=BothStrands())
    best = scan(model, batch; strands=BestStrand())
    lengths = Mimosa.scan_result_lengths(model, batch)

    @test nrows(scores) == 3
    @test rowlength(scores, 1) == expected_length
    @test rowlength(scores, 3) == 0
    @test nrows(both.forward) == 3
    @test nrows(both.reverse) == 3
    @test nrows(best) == 3
    @test lengths == [expected_length, npositions(model, length(sequences[2])), 0]
end

function _copy_motif_collection(directory, examples)
    collection = joinpath(directory, "motifs")
    mkpath(collection)
    for filename in ("foxa2.meme", "gata2.meme", "gata4.meme")
        cp(joinpath(examples, filename), joinpath(collection, filename))
    end
    return collection
end

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
