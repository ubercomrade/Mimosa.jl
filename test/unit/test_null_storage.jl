using Test
using Mimosa

@testset "Null distribution storage round-trip" begin
    # Create a null distribution
    gev = GEVFit(0.0, 0.5, 1.2, true, 42, -100.0)
    raw_scores = Float64[1.0, 2.0, 3.0, 4.0, 5.0, 1.5, 2.5, 3.5]
    dist = NullDistribution(
        "profile",
        "co",
        gev,
        raw_scores,
        [NullPair("m1", "t1", score) for score in raw_scores],
        8,
        2,
        [(query="m1", reason="too few")],
        "abc123",
        "def456",
        "seq_fp",
        "bg_fp",
    )

    mktempdir() do parent
        path = joinpath(parent, "null_bundle")
        savenull(path, dist)

        # Verify files exist
        @test isfile(joinpath(path, "manifest.toml"))
        @test isfile(joinpath(path, "data", "raw_null_scores.npy"))

        # Load and verify
        loaded = loadnull(path)

        @test loaded.strategy == "profile"
        @test loaded.metric == "co"
        @test loaded.n_null == 8
        @test loaded.n_queries == 2
        @test loaded.raw_scores == raw_scores
        @test loaded.fit isa GEVFit
        @test loaded.fit.shape ≈ 0.0 atol = 1e-10
        @test loaded.fit.location ≈ 0.5 atol = 1e-10
        @test loaded.fit.scale ≈ 1.2 atol = 1e-10
        @test loaded.fit.converged == true
        @test loaded.fit.iterations == 42
        @test loaded.model_collection_fingerprint == "abc123"
        @test loaded.relation_fingerprint == "def456"
        @test loaded.sequence_fingerprint == "seq_fp"
        @test loaded.background_fingerprint == "bg_fp"
        @test length(loaded.skipped) == 1
        @test loaded.skipped[1].query == "m1"
        @test loaded.skipped[1].reason == "too few"
    end
end

@testset "Null storage checksum validation" begin
    gev = GEVFit(-0.1, 0.0, 1.0, true, 10, -50.0)
    dist = NullDistribution(
        "profile",
        "co",
        gev,
        Float64[0.1, 0.2, 0.3, 0.4, 0.5],
        [NullPair("m", "t", score) for score in 0.1:0.1:0.5],
        5,
        1,
        [],
        nothing,
        nothing,
        "none",
        "none",
    )

    mktempdir() do parent
        path = joinpath(parent, "null_bundle")
        savenull(path, dist)

        # Corrupt the NPY file
        npy_path = joinpath(path, "data", "raw_null_scores.npy")
        open(npy_path, "a") do io
            write(io, UInt8(0xFF))
        end

        @test_throws ModelFormatError loadnull(path)
    end
end

@testset "Null storage format validation" begin
    mktempdir() do path
        # Missing manifest
        @test_throws ModelFormatError loadnull(path)

        # Wrong format
        mkpath(joinpath(path, "data"))
        open(joinpath(path, "manifest.toml"), "w") do io
            write(
                io, "format = \"other\"\nformat_version = 1\nkind = \"null_distribution\"\n"
            )
        end
        @test_throws ModelFormatError loadnull(path)
    end

    gev = GEVFit(0.0, 0.0, 1.0, true, 1, 0.0)
    dist = NullDistribution(
        "profile",
        "co",
        gev,
        Float64[1.0],
        [NullPair("m", "t", 1.0)],
        1,
        1,
        [],
        nothing,
        nothing,
        "none",
        "none",
    )
    mktempdir() do parent
        path = joinpath(parent, "legacy_null_bundle")
        savenull(path, dist)
        manifest_path = joinpath(path, "manifest.toml")
        manifest = read(manifest_path, String)
        write(
            manifest_path, replace(manifest, "format_version = 4" => "format_version = 3")
        )
        @test_throws ModelFormatError loadnull(path)
    end
end

@testset "Null storage: hostile manifest validation" begin
    gev = GEVFit(0.0, 0.5, 1.2, true, 42, -100.0)
    dist = NullDistribution(
        "profile",
        "co",
        gev,
        Float64[1.0, 2.0, 3.0],
        [NullPair("m", "t", score) for score in 1.0:1.0:3.0],
        3,
        1,
        NamedTuple{(:query, :reason),Tuple{String,String}}[],
        nothing,
        nothing,
        "none",
        "none",
    )

    mktempdir() do parent
        path = joinpath(parent, "null_bundle")
        savenull(path, dist)
        manifest_path = joinpath(path, "manifest.toml")
        original = read(manifest_path, String)
        checksum = match(r"checksum = \"(sha256:[0-9a-f]+)\"", original).captures[1]

        for bad_path in ["../outside.npy", "/tmp/outside.npy", raw"..\outside.npy"]
            write(manifest_path, replace(original, "data/raw_null_scores.npy" => bad_path))
            @test_throws ModelFormatError loadnull(path)
        end

        write(manifest_path, replace(original, checksum => "sha256:"))
        @test_throws ModelFormatError loadnull(path)

        write(manifest_path, replace(original, "n_null = 3" => "n_null = 300000000"))
        @test_throws ModelFormatError loadnull(path)

        write(manifest_path, original)
        rm(manifest_path)
        @test_throws ModelFormatError loadnull(path)
    end
end
