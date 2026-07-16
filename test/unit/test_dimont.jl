using Test
using Mimosa

const DIMONT_FIXTURES = joinpath(@__DIR__, "..", "fixtures")

@testset "Dimont model contract" begin
    _test_context_model_contract(Dimont)
end

@testset "Dimont parsing" begin
    # Test reading exampleD-model-1.xml (span=0, length=13)
    path = joinpath(DIMONT_FIXTURES, "exampleD-model-1.xml")
    @test isfile(path)
    m = read_dimont(path)
    @test m.name == "exampleD-model-1"
    @test m.span == 0
    @test m.motif_length == 13
    @test size(m.representation) == (5, 13)

    # Test reading stat_dimont-model-1.xml (span=3, length=5)
    path2 = joinpath(DIMONT_FIXTURES, "stat_dimont-model-1.xml")
    @test isfile(path2)
    m2 = read_dimont(path2)
    @test m2.span == 3
    @test m2.motif_length == 5
    @test size(m2.representation) == (625, 5)

    # Test reading PEAKS036274 (span=3, length=10)
    path3 = joinpath(DIMONT_FIXTURES, "PEAKS036274_FOXA1_P35582_MACS2-model-1.xml")
    @test isfile(path3)
    m3 = read_dimont(path3)
    @test m3.span == 3
    @test m3.motif_length == 10
    @test size(m3.representation) == (625, 10)

    # File not found
    @test_throws MimosaError read_dimont("nonexistent.xml")
end

@testset "Dimont scanning contract" begin
    model = read_dimont(joinpath(DIMONT_FIXTURES, "stat_dimont-model-1.xml"))
    sequence = UInt8[0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3]
    _test_scan_contract(model, sequence)
    _test_batch_scan_contract(model, (sequence, reverse(sequence)))
end

@testset "Dimont span=0 equivalence to order-0 BaMM scan" begin
    # A span=0 Dimont should produce the same scanning geometry as an order=0 BaMM:
    # kmer=1, context=0, window=motif_length, n_terms=motif_length
    path = joinpath(DIMONT_FIXTURES, "exampleD-model-1.xml")
    m = read_dimont(path)

    @test Mimosa.kmer(m) == 1
    @test Mimosa.context_length(m) == 0
    @test Mimosa.window_size(m) == m.motif_length
    @test Mimosa.scan_width(m) == m.motif_length

    # Verify scanning produces finite results
    seq = UInt8[0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3]
    fwd = scan(m, seq; strands=ForwardOnly())
    @test length(fwd) == 20 - 13 + 1  # 8
    @test all(isfinite, fwd)
end
