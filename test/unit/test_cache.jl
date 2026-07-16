using Test
using Mimosa
using SHA

@testset "Cache construction" begin
    # Cache does not create directories on construction
    mkpath(mktempdir())
    dir = joinpath(mktempdir(), "test_cache")
    cache = Cache(dir)
    @test cache.directory == dir
    @test cache.enabled == true
    @test !isdir(dir)  # not created yet

    # Disabled cache
    disabled = Cache(dir; enabled=false)
    @test disabled.enabled == false
end

@testset "Cache key contract" begin
    dir = mktempdir()
    cache = Cache(dir)
    @test_throws ArgumentError cache_set(cache, "human key", UInt8[1])
    cache_set(cache, "keep-key", UInt8[1, 2, 3])
    stage = joinpath(dir, ".mimosa-cache-stage-orphan")
    mkpath(stage)
    write(joinpath(stage, "partial"), UInt8[9])
    sentinel = joinpath(dir, "notes.txt")
    write(sentinel, "keep")
    @test clearcache(cache) == 1
    @test !ispath(stage)
    @test isfile(sentinel)
end

@testset "Cache set/get/has" begin
    dir = mktempdir()
    cache = Cache(dir)

    key = "test123"
    data = UInt8[1, 2, 3, 4, 5]

    # Initially absent
    @test !cache_has(cache, key)
    @test cache_get(cache, key) === nothing

    # Set data
    cache_set(cache, key, data; metadata=Dict("algorithm" => "test"))

    # Now present
    @test cache_has(cache, key)
    @test cache_get(cache, key) == data

    # Metadata
    meta = cache_get_meta(cache, key)
    @test meta !== nothing
    @test meta["algorithm"] == "test"
    @test startswith(meta["checksum"], "sha256:")
    @test meta["size"] == 5

    cache_set(
        cache,
        key,
        data;
        metadata=Dict("format_version" => 99, "checksum" => "sha256:bad", "size" => 0),
    )
    meta = cache_get_meta(cache, key)
    @test meta["format_version"] == Mimosa.CACHE_FORMAT_VERSION
    @test meta["size"] == length(data)
    @test cache_has(cache, key)
end

@testset "Cache disabled" begin
    # Use a non-existent subdirectory
    dir = joinpath(mktempdir(), "disabled_cache")
    cache = Cache(dir; enabled=false)

    @test !cache_has(cache, "key")
    @test cache_get(cache, "key") === nothing
    @test cache_get_meta(cache, "key") === nothing

    # cache_set is a no-op and should not create the directory
    cache_set(cache, "key", UInt8[1, 2, 3])
    @test !isdir(dir)
    @test !cache_has(cache, "key")
end

@testset "Cache clearcache" begin
    dir = mktempdir()
    cache = Cache(dir)

    # Write some entries
    cache_set(cache, "key1", UInt8[1, 2, 3])
    cache_set(cache, "key2", UInt8[4, 5, 6])
    @test cache_has(cache, "key1")
    @test cache_has(cache, "key2")

    # Clear single entry
    removed = clearcache(cache, "key1")
    @test removed == 1  # one committed cache entry directory
    @test !cache_has(cache, "key1")
    @test cache_has(cache, "key2")

    # Clear all
    count = clearcache(cache)
    @test count >= 1  # key2's committed entry
    @test !cache_has(cache, "key2")
end

@testset "Cache key containment" begin
    cache = Cache(mktempdir())
    for key in ("../escape", "/tmp/escape", "a\\\\b", "a\0b", "C:/escape", ".", "..")
        @test_throws ArgumentError cache_has(cache, key)
        @test_throws ArgumentError cache_get(cache, key)
        @test_throws ArgumentError cache_get_meta(cache, key)
        @test_throws ArgumentError cache_set(cache, key, UInt8[1])
        @test_throws ArgumentError clearcache(cache, key)
    end

    cache_set(cache, "keep", UInt8[1, 2])
    sentinel = joinpath(cache.directory, "sentinel.bin")
    write(sentinel, "keep")
    unrelated = joinpath(cache.directory, "unrelated")
    mkpath(unrelated)
    @test clearcache(cache, "unrelated") == 0
    @test isdir(unrelated)
    @test clearcache(cache) == 1
    @test isfile(sentinel)
    @test isdir(unrelated)

    if Sys.isunix()
        outside = tempname()
        write(outside, UInt8[9])
        symlink(outside, joinpath(cache.directory, "escape"))
        @test_throws ArgumentError cache_has(cache, "escape")
        @test islink(joinpath(cache.directory, "escape"))
        rm(outside; force=true)
    end
end

@testset "Cache corruption recovery" begin
    dir = mktempdir()
    cache = Cache(dir)

    key = "corrupt"
    data = UInt8[1, 2, 3, 4, 5]
    cache_set(cache, key, data)
    @test cache_has(cache, key)

    # Corrupt the data file
    data_path = joinpath(dir, key, "data.bin")
    write(data_path, UInt8[0, 0, 0, 0, 0])

    # Should be detected as corrupted (checksum mismatch), replaceable, and clearable.
    @test !cache_has(cache, key)
    @test cache_get(cache, key) === nothing
    cache_set(cache, key, data)
    @test cache_get(cache, key) == data

    write(data_path, UInt8[0, 0, 0, 0, 0])
    @test clearcache(cache, key) == 1
    @test !ispath(joinpath(dir, key))

    cache_set(cache, key, data)
    rm(joinpath(dir, key, "meta.toml"))
    @test clearcache(cache) == 1
    @test !ispath(joinpath(dir, key))
end

@testset "Concurrent cache commits retain a complete entry" begin
    dir = mktempdir()
    cache = Cache(dir)
    values = [fill(UInt8(i), 1024) for i in 1:16]
    tasks = [Threads.@spawn cache_set(cache, "shared", value) for value in values]
    fetch.(tasks)

    stored = cache_get(cache, "shared")
    @test stored in values
    @test cache_has(cache, "shared")
end

@testset "Cache key determinism" begin
    dir = mktempdir()
    cache = Cache(dir)

    # Same inputs → same key
    key1 = cache_key(cache, "pwm_scan", "model_fp_abc", "seq_fp_def", "forward")
    key2 = cache_key(cache, "pwm_scan", "model_fp_abc", "seq_fp_def", "forward")
    @test key1 == key2

    # Different inputs → different key
    key3 = cache_key(cache, "pwm_scan", "model_fp_abc", "seq_fp_def", "reverse")
    @test key1 != key3

    # Key length is 16 chars
    @test length(key1) == 16
end

@testset "Content fingerprint stability" begin
    # String fingerprints
    @test content_fingerprint("hello") == content_fingerprint("hello")
    @test content_fingerprint("hello") != content_fingerprint("world")

    # Byte array fingerprints
    @test content_fingerprint(UInt8[1, 2, 3]) == content_fingerprint(UInt8[1, 2, 3])
    @test content_fingerprint(UInt8[1, 2, 3]) != content_fingerprint(UInt8[1, 2, 4])

    # Array fingerprints
    arr1 = Float32[1.0 2.0; 3.0 4.0]
    arr2 = Float32[1.0 2.0; 3.0 4.0]
    @test content_fingerprint(arr1) == content_fingerprint(arr2)
    arr3 = Float32[1.0 2.0; 3.0 5.0]
    @test content_fingerprint(arr1) != content_fingerprint(arr3)
end

@testset "Model fingerprint" begin
    weights = Float32[
        0.5 -0.5 0.3
        -0.3 0.7 -0.2
        0.1 0.1 0.8
        -0.2 0.3 -0.1
        -0.3 -0.3 -0.3
    ]
    bg = (0.25f0, 0.25f0, 0.25f0, 0.25f0)
    pwm1 = PWM("model_a", weights, bg)
    pwm2 = PWM("model_a", weights, bg)
    pwm3 = PWM("model_b", weights, bg)

    # Same model → same fingerprint
    @test model_fingerprint(pwm1) == model_fingerprint(pwm2)

    # Different name → different fingerprint
    @test model_fingerprint(pwm1) != model_fingerprint(pwm3)

    # Different weights → different fingerprint
    weights2 = copy(weights)
    weights2[1, 1] = 0.6f0
    pwm4 = PWM("model_a", weights2, bg)
    @test model_fingerprint(pwm1) != model_fingerprint(pwm4)
end

@testset "ScoreProfile fingerprint" begin
    first = ScoreProfile("same", RaggedArray(Float32[1, 2], [1, 3]))
    second = ScoreProfile("same", RaggedArray(Float32[1, 3], [1, 3]))
    shifted = ScoreProfile("same", RaggedArray(Float32[1, 2], [1, 2, 3]))
    @test content_fingerprint(first) != content_fingerprint(second)
    @test content_fingerprint(first) != content_fingerprint(shifted)
end

@testset "Sequence fingerprint" begin
    seqs = [encode_sequence("ACGT"), encode_sequence("TTGG")]
    data = UInt8[]
    offsets = [1]
    for s in seqs
        append!(data, s)
        push!(offsets, length(data) + 1)
    end
    batch = EncodedSequenceBatch(data, offsets)

    # Same batch → same fingerprint
    @test sequence_fingerprint(batch) == sequence_fingerprint(batch)

    # Different batch → different fingerprint
    seqs2 = [encode_sequence("ACGT"), encode_sequence("TTGA")]
    data2 = UInt8[]
    offsets2 = [1]
    for s in seqs2
        append!(data2, s)
        push!(offsets2, length(data2) + 1)
    end
    batch2 = EncodedSequenceBatch(data2, offsets2)
    @test sequence_fingerprint(batch) != sequence_fingerprint(batch2)
end

@testset "Model collection fingerprint" begin
    weights = Float32[
        0.5 -0.5
        -0.3 0.7
        0.1 0.1
        -0.2 0.3
        -0.3 -0.3
    ]
    bg = (0.25f0, 0.25f0, 0.25f0, 0.25f0)
    m1 = PWM("m1", weights, bg)
    m2 = PWM("m2", weights, bg)

    # Same collection (any order) → same fingerprint
    @test model_collection_fingerprint([m1, m2]) == model_collection_fingerprint([m2, m1])

    # Different collection → different fingerprint
    @test model_collection_fingerprint([m1]) != model_collection_fingerprint([m1, m2])

    profile = ScoreProfile("profile", build_ragged([Float32[0.1, 0.2]]))
    @test model_collection_fingerprint([profile]) ==
        model_collection_fingerprint(AbstractProfileSource[profile])
    @test model_collection_fingerprint([m1, profile]) ==
        model_collection_fingerprint([profile, m1])
end
