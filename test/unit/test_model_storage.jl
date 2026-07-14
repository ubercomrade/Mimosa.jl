using Test
using Mimosa

@testset "Model storage round-trip: PWM" begin
    weights = Float32[
        0.5 -0.5 0.3 0.1
        -0.3 0.7 -0.2 0.2
        0.1 0.1 0.8 -0.3
        -0.2 0.3 -0.1 0.4
        -0.3 -0.3 -0.3 -0.3
    ]
    bg = (0.25f0, 0.25f0, 0.30f0, 0.20f0)
    pwm = PWM("test_pwm", weights, bg)

    dir = mktempdir()
    bundle_path = joinpath(dir, "pwm_bundle")

    writemodel(bundle_path, pwm)
    loaded = readmodel(bundle_path)

    @test loaded isa PWM
    @test loaded.name == pwm.name
    @test loaded.weights == pwm.weights
    @test loaded.background == pwm.background
end

@testset "Model storage rejects legacy PFM bundles" begin
    dir = mktempdir()
    bundle_path = joinpath(dir, "pfm_bundle")
    mkpath(bundle_path)
    write(
        joinpath(bundle_path, "manifest.toml"),
        "format = \"mimosa\"\nformat_version = 1\nkind = \"pfm\"\nname = \"legacy\"\n",
    )
    @test_throws ModelFormatError readmodel(bundle_path)
end

@testset "Model storage round-trip: BaMM" begin
    # Order 1, 4 positions: 5^2 = 25 rows
    rep = Float32.(reshape(1.0:75, 25, 3))
    model = BaMM("test_bamm", rep, 1, 3)

    dir = mktempdir()
    bundle_path = joinpath(dir, "bamm_bundle")

    writemodel(bundle_path, model)
    loaded = readmodel(bundle_path)

    @test loaded isa BaMM
    @test loaded.name == model.name
    @test loaded.representation == model.representation
    @test loaded.order == model.order
    @test loaded.motif_length == model.motif_length
end

@testset "Model storage round-trip: SiteGA" begin
    rep = Float32.(reshape(1.0:75, 25, 3))
    model = SiteGA("test_sitega", rep, 3)

    dir = mktempdir()
    bundle_path = joinpath(dir, "sitega_bundle")

    writemodel(bundle_path, model)
    loaded = readmodel(bundle_path)

    @test loaded isa SiteGA
    @test loaded.name == model.name
    @test loaded.representation == model.representation
    @test loaded.motif_length == model.motif_length
end

@testset "Model storage round-trip: Dimont" begin
    # Span 1, 3 positions: 5^2 = 25 rows
    rep = Float32.(reshape(1.0:75, 25, 3))
    model = Dimont("test_dimont", rep, 1, 3)

    dir = mktempdir()
    bundle_path = joinpath(dir, "dimont_bundle")

    writemodel(bundle_path, model)
    loaded = readmodel(bundle_path)

    @test loaded isa Dimont
    @test loaded.name == model.name
    @test loaded.representation == model.representation
    @test loaded.span == model.span
    @test loaded.motif_length == model.motif_length
end

@testset "Model storage round-trip: Slim" begin
    rep = Float32.(reshape(1.0:75, 25, 3))
    model = Slim("test_slim", rep, 1, 3)

    dir = mktempdir()
    bundle_path = joinpath(dir, "slim_bundle")

    writemodel(bundle_path, model)
    loaded = readmodel(bundle_path)

    @test loaded isa Slim
    @test loaded.name == model.name
    @test loaded.representation == model.representation
    @test loaded.span == model.span
    @test loaded.motif_length == model.motif_length
end

@testset "Model storage: checksum validation" begin
    weights = Float32[
        0.5 -0.5
        -0.3 0.7
        0.1 0.1
        -0.2 0.3
        -0.3 -0.3
    ]
    bg = (0.25f0, 0.25f0, 0.25f0, 0.25f0)
    pwm = PWM("test", weights, bg)

    dir = mktempdir()
    bundle_path = joinpath(dir, "pwm_bundle")

    writemodel(bundle_path, pwm)

    # Corrupt the data file
    data_file = joinpath(bundle_path, "data", "weights.bin")
    write(data_file, UInt8[0, 0, 0, 0, 0])

    @test_throws MimosaError readmodel(bundle_path)
end

@testset "Model storage: unknown format" begin
    dir = mktempdir()
    bundle_path = joinpath(dir, "bad_bundle")
    mkpath(joinpath(bundle_path, "data"))

    # Write manifest with wrong format
    manifest_path = joinpath(bundle_path, "manifest.toml")
    write(
        manifest_path,
        """
format = "wrong"
format_version = 2
kind = "pwm"
name = "test"
""",
    )

    @test_throws MimosaError readmodel(bundle_path)
end

@testset "Model storage: version too high" begin
    weights = Float32[
        0.5 -0.5
        -0.3 0.7
        0.1 0.1
        -0.2 0.3
        -0.3 -0.3
    ]
    bg = (0.25f0, 0.25f0, 0.25f0, 0.25f0)
    pwm = PWM("test", weights, bg)

    dir = mktempdir()
    bundle_path = joinpath(dir, "pwm_bundle")

    writemodel(bundle_path, pwm)

    # Bump the version in manifest
    manifest_path = joinpath(bundle_path, "manifest.toml")
    content = read(manifest_path, String)
    content = replace(content, "format_version = 2" => "format_version = 999")
    write(manifest_path, content)

    @test_throws MimosaError readmodel(bundle_path)
end

@testset "readmodel: legacy fallback still works" begin
    # Ensure legacy file reading still works alongside bundle reading
    # by reading a .pfm file (non-directory path)
    dir = mktempdir()
    pfm_path = joinpath(dir, "test.pfm")
    write(
        pfm_path,
        """
0.3 0.2 0.4 0.1
0.2 0.3 0.1 0.4
0.1 0.4 0.2 0.3
0.4 0.1 0.3 0.2
""",
    )

    model = readmodel(pfm_path)
    @test model isa PWM
    @test size(model.weights) == (5, 4)
end

function _storage_test_pwm()
    weights = Float32[
        0.5 -0.5
        -0.3 0.7
        0.1 0.1
        -0.2 0.3
        -0.3 -0.3
    ]
    return PWM("hostile_test", weights, (0.25f0, 0.25f0, 0.25f0, 0.25f0))
end

function _refresh_model_checksum!(manifest_path, data_path)
    content = read(manifest_path, String)
    old = match(r"checksum = \"(sha256:[0-9a-f]+)\"", content)
    old === nothing && error("test fixture has no checksum")
    updated = replace(
        content, old.captures[1] => "sha256:" * Mimosa._file_sha256(data_path)
    )
    write(manifest_path, updated)
    return nothing
end

@testset "Model storage: hostile bundle validation" begin
    mktempdir() do root
        bundle = joinpath(root, "model_bundle")
        writemodel(bundle, _storage_test_pwm())
        manifest_path = joinpath(bundle, "manifest.toml")
        data_path = joinpath(bundle, "data", "weights.bin")
        original_manifest = read(manifest_path, String)
        original_data = read(data_path)

        for bad_path in [
            "../outside.bin",
            "/tmp/outside.bin",
            raw"..\outside.bin",
            raw"C:\tmp\outside.bin",
            raw"C:/tmp/outside.bin",
        ]
            write(manifest_path, replace(original_manifest, "data/weights.bin" => bad_path))
            @test_throws ModelFormatError readmodel(bundle)
        end

        checksum = match(r"checksum = \"(sha256:[0-9a-f]+)\"", original_manifest).captures[1]
        for bad_checksum in [
            "sha256:",
            "sha256:" * repeat("A", 64),
            "md5:" * repeat("0", 32),
            "sha256:" * repeat("0", 63),
        ]
            write(manifest_path, replace(original_manifest, checksum => bad_checksum))
            @test_throws ModelFormatError readmodel(bundle)
        end

        write(manifest_path, replace(original_manifest, "checksum = \"$checksum\"\n" => ""))
        @test_throws ModelFormatError readmodel(bundle)

        for version in ["0", "-1", "1", "1.0", "\"2\""]
            write(
                manifest_path,
                replace(
                    original_manifest, "format_version = 2" => "format_version = $version"
                ),
            )
            @test_throws ModelFormatError readmodel(bundle)
        end

        write(
            manifest_path,
            replace(original_manifest, "shape = [5, 2]" => "shape = [100000001, 2]"),
        )
        @test_throws ModelFormatError readmodel(bundle)

        write(
            manifest_path,
            replace(original_manifest, "name = \"hostile_test\"" => "name = [1]"),
        )
        @test_throws ModelFormatError readmodel(bundle)

        write(manifest_path, "format = [\n")
        @test_throws ModelFormatError readmodel(bundle)

        write(manifest_path, repeat("# oversized manifest\n", 60_000))
        @test_throws ModelFormatError readmodel(bundle)

        write(manifest_path, original_manifest)
        bad_data = copy(original_data)
        bad_data[1:4] .= 0xff
        write(data_path, bad_data)
        _refresh_model_checksum!(manifest_path, data_path)
        @test_throws ModelFormatError readmodel(bundle)

        write(data_path, original_data[1:(end - 1)])
        _refresh_model_checksum!(manifest_path, data_path)
        @test_throws ModelFormatError readmodel(bundle)

        write(data_path, original_data)
        _refresh_model_checksum!(manifest_path, data_path)
        @test readmodel(bundle) isa PWM

        write(data_path, vcat(original_data, UInt8[0xff]))
        _refresh_model_checksum!(manifest_path, data_path)
        @test_throws ModelFormatError readmodel(bundle)

        write(data_path, original_data)
        _refresh_model_checksum!(manifest_path, data_path)
        if !Sys.iswindows()
            outside = joinpath(root, "outside.bin")
            write(outside, original_data)
            rm(data_path)
            symlink(outside, data_path)
            @test_throws ModelFormatError readmodel(bundle)
            rm(data_path)
            rm(bundle; recursive=true)
            writemodel(bundle, _storage_test_pwm())
        end

        @test isempty(
            filter(name -> startswith(name, ".model_bundle.mimosa-stage-"), readdir(root))
        )
    end
end

@testset "Model storage: staged writes leave complete bundles" begin
    mktempdir() do root
        bundle = joinpath(root, "model_bundle")
        writemodel(bundle, _storage_test_pwm())
        @test readmodel(bundle) isa PWM
        write(joinpath(root, "not_a_directory"), "occupied")
        @test_throws InvariantError writemodel(
            joinpath(root, "not_a_directory"), _storage_test_pwm()
        )
        @test isempty(
            filter(
                name -> startswith(name, ".not_a_directory.mimosa-stage-"), readdir(root)
            ),
        )
    end
end
