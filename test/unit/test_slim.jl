using Test
using Mimosa

const SLIM_FIXTURES = joinpath(@__DIR__, "..", "fixtures", "slim")

@testset "Slim model contract" begin
    _test_context_model_contract(Slim)
end

@testset "Slim parsing" begin
    # Test reading example-model-1.xml (span=5, length=15)
    path = joinpath(SLIM_FIXTURES, "example-model-1.xml")
    @test isfile(path)
    m = read_slim(path)
    @test m.name == "example-model-1"
    @test m.span == 5
    @test m.motif_length == 15
    @test size(m.representation) == (15625, 15)

    # Test reading PEAKS036274 (span=5, length=12)
    path2 = joinpath(SLIM_FIXTURES, "PEAKS036274_FOXA1_P35582_MACS2-model-2.xml")
    @test isfile(path2)
    m2 = read_slim(path2)
    @test m2.span == 5
    @test m2.motif_length == 12
    @test size(m2.representation) == (15625, 12)

    # File not found
    @test_throws MimosaError read_slim("nonexistent.xml")
end

@testset "Slim scanning contract" begin
    model = read_slim(joinpath(SLIM_FIXTURES, "example-model-1.xml"))
    sequence = UInt8[
        0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0
    ]
    _test_scan_contract(model, sequence)
    _test_batch_scan_contract(model, (sequence, reverse(sequence)))
end

@testset "Slim span=0 equivalence to order-0 BaMM scan geometry" begin
    # A span=0 Slim should produce the same scanning geometry as an order=0 BaMM:
    # kmer=1, context=0, window=motif_length, n_terms=motif_length
    rep = Matrix{Float32}(undef, 5, 8)
    fill!(rep, 0.0f0)
    m = Slim("span0", rep, 0, 8)

    @test Mimosa.kmer(m) == 1
    @test Mimosa.context_length(m) == 0
    @test Mimosa.window_size(m) == m.motif_length
    @test Mimosa.scan_width(m) == m.motif_length
end
