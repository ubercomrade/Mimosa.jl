using Test

const REPO_ROOT = dirname(dirname(@__DIR__))
const MIMOSA_APP = joinpath(REPO_ROOT, "app", "mimosa.jl")
const EXAMPLES = joinpath(REPO_ROOT, "examples")

# Helper: run a command and return the exit code, suppressing pipeline errors
# for non-zero exit codes (which are expected in error-case tests).
function _run_exitcode(cmd::Cmd)
    try
        run(pipeline(cmd; stderr=devnull, stdout=devnull); wait=true)
        return 0
    catch e
        if e isa ProcessFailedException
            return first(e.procs).exitcode
        else
            rethrow(e)
        end
    end
end

@testset "CLI subprocess: help and version" begin
    code = _run_exitcode(`$(Base.julia_cmd()) --project=$(REPO_ROOT) $(MIMOSA_APP) --help`)
    @test code == 0

    code = _run_exitcode(
        `$(Base.julia_cmd()) --project=$(REPO_ROOT) $(MIMOSA_APP) --version`
    )
    @test code == 0

    code = _run_exitcode(`$(Base.julia_cmd()) --project=$(REPO_ROOT) $(MIMOSA_APP) -h`)
    @test code == 0

    code = _run_exitcode(`$(Base.julia_cmd()) --project=$(REPO_ROOT) $(MIMOSA_APP) -V`)
    @test code == 0

    code = _run_exitcode(
        `$(Base.julia_cmd()) --project=$(REPO_ROOT) -m Mimosa.CLIApp --version`
    )
    @test code == 0
end

@testset "CLI subprocess: no args exits 1" begin
    code = _run_exitcode(`$(Base.julia_cmd()) --project=$(REPO_ROOT) $(MIMOSA_APP)`)
    @test code == 1
end

@testset "CLI subprocess: profile comparison" begin
    query = joinpath(EXAMPLES, "pif4.meme")
    target = joinpath(EXAMPLES, "gata2.meme")
    code = _run_exitcode(
        `$(Base.julia_cmd()) --project=$(REPO_ROOT) $(MIMOSA_APP) profile $(query) $(target) --model1-type pwm --model2-type pwm --metric co --num-sequences 50 --seq-length 100 --seed 42`,
    )
    @test code == 0
end

@testset "CLI subprocess: cache clear" begin
    dir = mktempdir()
    cache_dir = joinpath(dir, "mimosa-cache")
    code = _run_exitcode(
        `$(Base.julia_cmd()) --project=$(REPO_ROOT) $(MIMOSA_APP) cache clear --cache-dir $(cache_dir)`,
    )
    @test code == 0
end

@testset "CLI subprocess: build-null profile strategy" begin
    dir = mktempdir()
    coll_dir = _copy_motif_collection(dir, EXAMPLES)
    groups_path = joinpath(dir, "groups.tsv")
    write(groups_path, "motif\tgroup\nMA0047.3\tA\nMA0036.2\tB\nMA0482.2\tC\n")
    output_path = joinpath(dir, "null")

    code = _run_exitcode(
        `$(Base.julia_cmd()) --project=$(REPO_ROOT) $(MIMOSA_APP) build-null $(coll_dir) --model-type pwm --groups $(groups_path) --metric co --num-sequences 50 --seq-length 100 --seed 42 --output $(output_path)`,
    )
    @test code == 0
    @test isfile(joinpath(output_path, "manifest.toml"))
end

@testset "CLI subprocess: unknown command exits 1" begin
    code = _run_exitcode(
        `$(Base.julia_cmd()) --project=$(REPO_ROOT) $(MIMOSA_APP) bogus-command`
    )
    @test code == 1
end

@testset "CLI subprocess: missing required file exits 2" begin
    code = _run_exitcode(
        `$(Base.julia_cmd()) --project=$(REPO_ROOT) $(MIMOSA_APP) profile /nonexistent.meme /also-nonexistent.meme --model1-type pwm --model2-type pwm`,
    )
    @test code == 2
end
