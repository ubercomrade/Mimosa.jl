using Test
using Mimosa

function _stored_null(raw_scores=Float64[1.0, 2.0, 3.0, 4.0, 5.0])
    gev = GEVFit(0.0, 0.5, 1.2, true, 42, -100.0)
    return NullDistribution(
        "profile",
        "co",
        gev,
        raw_scores,
        [NullPair("m1", "m2", score) for score in raw_scores],
        length(raw_scores),
        3,
        "pwm",
        true,
        42,
        "random-ordered-pairs-v1",
        "abc123",
        "seq_fp",
        "bg_fp",
    )
end

@testset "Null distribution storage round-trip" begin
    raw_scores = Float64[1.0, 2.0, 3.0, 4.0, 5.0, 1.5, 2.5, 3.5]
    dist = _stored_null(raw_scores)

    mktempdir() do parent
        path = joinpath(parent, "null_bundle")
        savenull(path, dist)

        @test isfile(joinpath(path, "manifest.toml"))
        @test isfile(joinpath(path, "data", "raw_null_scores.npy"))

        loaded = loadnull(path)
        @test loaded.strategy == "profile"
        @test loaded.metric == "co"
        @test loaded.n_null == 8
        @test loaded.n_models == 3
        @test loaded.model_type == "pwm"
        @test loaded.shuffle
        @test loaded.seed == 42
        @test loaded.sampling_version == "random-ordered-pairs-v1"
        @test loaded.raw_scores == raw_scores
        @test loaded.fit isa GEVFit
        @test loaded.fit.shape ≈ 0.0 atol = 1e-10
        @test loaded.fit.location ≈ 0.5 atol = 1e-10
        @test loaded.fit.scale ≈ 1.2 atol = 1e-10
        @test loaded.fit.converged
        @test loaded.fit.iterations == 42
        @test loaded.model_collection_fingerprint == "abc123"
        @test loaded.sequence_fingerprint == "seq_fp"
        @test loaded.background_fingerprint == "bg_fp"
    end
end

@testset "Null storage checksum validation" begin
    dist = _stored_null()
    mktempdir() do parent
        path = joinpath(parent, "null_bundle")
        savenull(path, dist)
        open(joinpath(path, "data", "raw_null_scores.npy"), "a") do io
            write(io, UInt8(0xFF))
        end
        @test_throws ModelFormatError loadnull(path)
    end
end

@testset "Null storage format validation" begin
    mktempdir() do path
        @test_throws ModelFormatError loadnull(path)
        mkpath(joinpath(path, "data"))
        open(joinpath(path, "manifest.toml"), "w") do io
            write(
                io, "format = \"other\"\nformat_version = 1\nkind = \"null_distribution\"\n"
            )
        end
        @test_throws ModelFormatError loadnull(path)
    end

    mktempdir() do parent
        path = joinpath(parent, "legacy_null_bundle")
        savenull(path, _stored_null())
        manifest_path = joinpath(path, "manifest.toml")
        manifest = read(manifest_path, String)
        write(
            manifest_path, replace(manifest, "format_version = 5" => "format_version = 4")
        )
        @test_throws ModelFormatError loadnull(path)
    end
end

@testset "Null storage: hostile manifest validation" begin
    dist = _stored_null(Float64[1.0, 2.0, 3.0])
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
