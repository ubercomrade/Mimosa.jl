using Test
using Mimosa

const REPO_ROOT = dirname(dirname(@__DIR__))
const EXAMPLES = joinpath(REPO_ROOT, "examples")

@testset "CLI main: help and version" begin
    # Global help
    @test Mimosa.main(["--help"]) == 0
    @test Mimosa.main(["--version"]) == 0
    @test Mimosa.main(["-h"]) == 0
    @test Mimosa.main(["-V"]) == 0

    # No args → help to stderr, exit 1
    @test Mimosa.main(String[]) == 1
end

@testset "CLI main: unknown command" begin
    @test Mimosa.main(["bogus"]) == 1
end

@testset "CLI rejects removed commands" begin
    for command in ("motif", "inspect-model", "convert-model")
        @test Mimosa.main([command, "--help"]) == 1
    end
end

@testset "CLI profile: score profile comparison" begin
    query = joinpath(EXAMPLES, "scores_1.fasta")
    target = joinpath(EXAMPLES, "scores_2.fasta")
    for metric in ("co", "co_rowwise", "dice", "dice_rowwise", "cosine")
        code = Mimosa.main([
            "profile",
            query,
            target,
            "--model1-type",
            "scores",
            "--model2-type",
            "scores",
            "--metric",
            metric,
        ])
        @test code == 0
    end
end

@testset "CLI profile: help" begin
    @test Mimosa.main(["profile", "--help"]) == 0
end

@testset "CLI profile: missing required args" begin
    code = Mimosa.main(["profile", "a.fasta", "b.fasta"])
    @test code == 1

    code = Mimosa.main(["profile", "--model1-type", "scores", "--model2-type", "scores"])
    @test code == 1
end

@testset "CLI profile: motif-derived with random sequences" begin
    query = joinpath(EXAMPLES, "pif4.meme")
    target = joinpath(EXAMPLES, "pif4.meme")
    code = Mimosa.main([
        "profile",
        query,
        target,
        "--model1-type",
        "pwm",
        "--model2-type",
        "pwm",
        "--metric",
        "co",
        "--num-sequences",
        "50",
        "--seq-length",
        "100",
        "--seed",
        "42",
    ])
    @test code == 0
end

@testset "CLI profile: motif-derived with FASTA" begin
    query = joinpath(EXAMPLES, "pif4.meme")
    fasta = joinpath(EXAMPLES, "foreground.fa")
    code = Mimosa.main([
        "profile",
        query,
        query,
        "--model1-type",
        "pwm",
        "--model2-type",
        "pwm",
        "--metric",
        "co",
        "--fasta",
        fasta,
    ])
    @test code == 0
end

@testset "CLI profile: BaMM model-derived" begin
    query = joinpath(EXAMPLES, "foxa2.ihbcp")
    target = joinpath(EXAMPLES, "gata2.ihbcp")
    code = Mimosa.main([
        "profile",
        query,
        target,
        "--model1-type",
        "bamm",
        "--model2-type",
        "bamm",
        "--metric",
        "co",
        "--num-sequences",
        "50",
        "--seq-length",
        "100",
        "--seed",
        "42",
    ])
    @test code == 0
end

@testset "CLI cache: clear empty cache" begin
    dir = mktempdir()
    cache_dir = joinpath(dir, "mimosa-cache")
    code = Mimosa.main(["cache", "clear", "--cache-dir", cache_dir])
    @test code == 0
end

@testset "CLI cache: clear with entries" begin
    dir = mktempdir()
    cache_dir = joinpath(dir, "mimosa-cache")
    cache = Cache(cache_dir)
    cache_set(cache, "test-key", b"test-data")
    @test cache_has(cache, "test-key")

    code = Mimosa.main(["cache", "clear", "--cache-dir", cache_dir])
    @test code == 0
    @test !cache_has(cache, "test-key")
end

@testset "CLI cache: help" begin
    @test Mimosa.main(["cache", "--help"]) == 0
end

@testset "CLI cache: invalid subcommand" begin
    code = Mimosa.main(["cache", "bogus"])
    @test code == 1
end

@testset "CLI build-null: help" begin
    @test Mimosa.main(["build-null", "--help"]) == 0
end

@testset "CLI build-null: missing required args" begin
    code = Mimosa.main(["build-null"])
    @test code == 1

    code = Mimosa.main(["build-null", "motifs.meme"])
    @test code == 1
end

@testset "CLI build-null: profile strategy with random sequences" begin
    dir = mktempdir()
    # Use multiple MEME files as motif collection (directory)
    coll_dir = joinpath(dir, "motifs")
    mkpath(coll_dir)
    cp(joinpath(EXAMPLES, "foxa2.meme"), joinpath(coll_dir, "foxa2.meme"))
    cp(joinpath(EXAMPLES, "gata2.meme"), joinpath(coll_dir, "gata2.meme"))
    cp(joinpath(EXAMPLES, "gata4.meme"), joinpath(coll_dir, "gata4.meme"))

    # Create groups TSV with motif names in different groups
    groups_path = joinpath(dir, "groups.tsv")
    write(groups_path, "motif\tgroup\nMA0047.3\tA\nMA0036.2\tB\nMA0482.2\tC\n")

    output_path = joinpath(dir, "null")

    code = Mimosa.main([
        "build-null",
        coll_dir,
        "--model-type",
        "pwm",
        "--groups",
        groups_path,
        "--metric",
        "co",
        "--num-sequences",
        "100",
        "--seq-length",
        "100",
        "--output",
        output_path,
    ])
    @test code == 0
    # Null output is a directory with manifest.toml
    @test isfile(joinpath(output_path, "manifest.toml"))

    # Annotation accepts only a bundle built for the executed strategy and metric.
    code = Mimosa.main([
        "profile",
        joinpath(EXAMPLES, "foxa2.meme"),
        joinpath(EXAMPLES, "gata2.meme"),
        "--model1-type",
        "pwm",
        "--model2-type",
        "pwm",
        "--metric",
        "co",
        "--num-sequences",
        "100",
        "--seq-length",
        "100",
        "--pvalue",
        "--null-distribution",
        output_path,
        "--effective-number-of-targets",
        "3",
    ])
    @test code == 0

    code = Mimosa.main([
        "profile",
        joinpath(EXAMPLES, "foxa2.meme"),
        joinpath(EXAMPLES, "gata2.meme"),
        "--model1-type",
        "pwm",
        "--model2-type",
        "pwm",
        "--metric",
        "dice",
        "--num-sequences",
        "100",
        "--seq-length",
        "100",
        "--pvalue",
        "--null-distribution",
        output_path,
    ])
    @test code == 1
end

@testset "CLI build-null: profile strategy passes FASTA and metric" begin
    dir = mktempdir()
    coll_dir = _copy_motif_collection(dir, EXAMPLES)
    groups_path = joinpath(dir, "groups.tsv")
    write(groups_path, "motif\tgroup\nMA0047.3\tA\nMA0036.2\tB\nMA0482.2\tC\n")
    output_path = joinpath(dir, "profile_null")

    code = Mimosa.main([
        "build-null",
        coll_dir,
        "--model-type",
        "pwm",
        "--groups",
        groups_path,
        "--metric",
        "dice",
        "--fasta",
        joinpath(EXAMPLES, "foreground.fa"),
        "--search-range",
        "2",
        "--window-radius",
        "2",
        "--realign-window",
        "1",
        "--min-logfpr",
        "-2.0",
        "--output",
        output_path,
    ])
    @test code == 0
    dist = loadnull(output_path)
    @test dist.strategy == "profile"
    @test dist.metric == "dice"
    @test dist.sequence_fingerprint != "none"

    code = Mimosa.main([
        "profile",
        joinpath(EXAMPLES, "foxa2.meme"),
        joinpath(EXAMPLES, "gata2.meme"),
        "--model1-type",
        "pwm",
        "--model2-type",
        "pwm",
        "--metric",
        "dice",
        "--fasta",
        joinpath(EXAMPLES, "foreground.fa"),
        "--search-range",
        "2",
        "--window-radius",
        "2",
        "--realign-window",
        "1",
        "--min-logfpr",
        "-2.0",
        "--pvalue",
        "--null-distribution",
        output_path,
    ])
    @test code == 0
end

@testset "CLI build-null: removed strategy option" begin
    code = Mimosa.main([
        "build-null",
        "motifs.meme",
        "--model-type",
        "pwm",
        "--groups",
        "groups.tsv",
        "--strategy",
        "bogus",
        "--output",
        "out",
    ])
    @test code == 1
end

@testset "CLI build-null: invalid profile metric" begin
    code = Mimosa.main([
        "build-null",
        "motifs.meme",
        "--model-type",
        "pwm",
        "--groups",
        "groups.tsv",
        "--metric",
        "pcc",
        "--output",
        "out",
    ])
    @test code == 1
end

@testset "CLI build-null: threads alias works" begin
    dir = mktempdir()
    coll_dir = _copy_motif_collection(dir, EXAMPLES)
    groups_path = joinpath(dir, "groups.tsv")
    write(groups_path, "motif\tgroup\nMA0047.3\tA\nMA0036.2\tB\nMA0482.2\tC\n")
    output_path = joinpath(dir, "null_threaded")

    # --jobs alias should work same as --threads
    code = Mimosa.main([
        "build-null",
        coll_dir,
        "--model-type",
        "pwm",
        "--groups",
        groups_path,
        "--metric",
        "co",
        "--output",
        output_path,
        "--jobs",
        "1",
    ])
    @test code == 0
    @test isfile(joinpath(output_path, "manifest.toml"))
end

@testset "Parallel null build: serial == threaded" begin
    weights1 = Float32[
        0.5 -0.5 0.3
        -0.3 0.7 -0.2
        0.1 0.1 0.8
        -0.2 0.3 -0.1
        -0.3 -0.3 -0.3
    ]
    weights2 = Float32[
        0.3 0.2 0.5
        0.1 0.8 0.1
        0.2 0.3 0.4
        0.1 0.1 0.2
        -0.1 -0.1 -0.1
    ]
    weights3 = Float32[
        0.4 -0.4 0.2
        -0.2 0.6 -0.1
        0.2 0.2 0.7
        -0.1 0.4 0.0
        -0.2 -0.2 -0.2
    ]
    bg = (0.25f0, 0.25f0, 0.25f0, 0.25f0)
    m1 = PWM("m1", weights1, bg)
    m2 = PWM("m2", weights2, bg)
    m3 = PWM("m3", weights3, bg)

    models = [m1, m2, m3]

    tsv_content = "motif\tgroup\nm1\tA\nm2\tB\nm3\tC\n"
    dir = mktempdir()
    rel_path = joinpath(dir, "relations.tsv")
    write(rel_path, tsv_content)
    relations = parse_group_relations(rel_path; known_names=Set(["m1", "m2", "m3"]))

    sequences = make_random_sequences(2, 30; seed=7)
    serial_result = build_null(
        models, relations; sequences=sequences, execution=SerialExecution()
    )
    serial_scores = serial_result.distribution.raw_scores

    for nt in (1, 2, 4)
        threaded_result = build_null(
            models, relations; sequences=sequences, execution=ThreadedExecution(nt)
        )
        threaded_scores = threaded_result.distribution.raw_scores
        @test threaded_scores == serial_scores
        @test threaded_result.total_comparisons == serial_result.total_comparisons
    end
end

@testset "Model storage: cross-format compatibility" begin
    pwm = readmodel(joinpath(EXAMPLES, "pif4.meme"))
    dir = mktempdir()
    bundle = joinpath(dir, "pif4_bundle")
    writemodel(bundle, pwm)
    loaded = readmodel(bundle)

    seq = encode_sequence("ACGTACGTACGTACGTACGTACGTAC")
    original_scan = scan(pwm, seq; strands=ForwardOnly())
    loaded_scan = scan(loaded, seq; strands=ForwardOnly())
    @test original_scan == loaded_scan
end
