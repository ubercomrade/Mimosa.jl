using Test
using Mimosa
using SHA

@testset "Execution" begin
    @test Execution() isa Execution
    @test Execution() == Execution(1)
    @test Execution(4).ntasks == 4
    @test Execution(1).ntasks == 1
    @test_throws ArgumentError Execution(0)

    # show
    @test sprint(show, Execution()) == "Execution(ntasks=1)"
    @test sprint(show, Execution(2)) == "Execution(ntasks=2)"
end

@testset "CLI execution" begin
    parsed = Mimosa.CLIParsed("profile")
    parsed.options["threads"] = string(Threads.nthreads())
    execution = Mimosa._execution(parsed)
    @test execution == Execution(Threads.nthreads())

    parsed.options["threads"] = string(Threads.nthreads() + 1)
    @test_throws Mimosa.CLIError Mimosa._execution(parsed)
end

@testset "CLI progress renderer" begin
    io = IOBuffer()
    progress = ProgressBar(io; width=5, refresh_seconds=0)
    progress((; stage=:prepare, current=0, total=2, label=""))
    progress((; stage=:prepare, current=1, total=2, label="m1"))
    progress((; stage=:prepare, current=2, total=2, label="m2"))
    progress((; stage=:null, current=0, total=3, label=""))
    progress((; stage=:null, current=3, total=3, label="m3"))
    output = String(take!(io))
    @test occursin("Preparing models", output)
    @test occursin("Building null", output)
    @test occursin("2/2", output)
    @test occursin("3/3", output)
    @test !progress.active

    parsed = Mimosa.CLIParsed("build-null")
    push!(parsed.flags, "quiet")
    @test Mimosa._cli_progress(parsed) === nothing
    @test !Mimosa._is_terminal_output(IOBuffer())
end

@testset "_parallel_chunks" begin
    for execution in (Execution(), Execution(4))
        visits = zeros(Int, 37)
        chunk_ids = zeros(Int, 37)
        Mimosa._parallel_chunks(execution, length(visits)) do first, last, chunk
            for i in first:last
                visits[i] += 1
                chunk_ids[i] = chunk
            end
        end
        @test visits == ones(Int, length(visits))
        @test issorted(chunk_ids)
        @test maximum(chunk_ids) <= min(execution.ntasks, Threads.nthreads())
    end

    Mimosa._parallel_chunks(Execution(4), 0) do _, _, _
        @test false
    end
end

@testset "_parallel_for serial" begin
    # Basic serial execution
    results = Vector{Int}(undef, 10)
    Mimosa._parallel_for(Execution(), 10) do i
        results[i] = i * 2
    end
    @test results == collect(2:2:20)
end

@testset "_parallel_for threaded" begin
    # Threaded execution should produce same results as serial
    for ntasks in (1, 2, 4)
        results_t = Vector{Int}(undef, 100)
        Mimosa._parallel_for(Execution(ntasks), 100) do i
            results_t[i] = i^2
        end
        results_s = Vector{Int}(undef, 100)
        Mimosa._parallel_for(Execution(), 100) do i
            results_s[i] = i^2
        end
        @test results_t == results_s
    end

    # Edge case: n=1
    results = Vector{Int}(undef, 1)
    Mimosa._parallel_for(Execution(4), 1) do i
        results[i] = 42
    end
    @test results == [42]

    # Edge case: n=0
    Mimosa._parallel_for(Execution(4), 0) do i
        @test false  # should not execute
    end
end

@testset "bounded and weighted execution" begin
    @test Mimosa._effective_ntasks(Execution(100), 100) == Threads.nthreads()

    visits = zeros(Int, 37)
    costs = [i % 7 == 0 ? 1000 : 1 for i in eachindex(visits)]
    Mimosa._parallel_for_weighted(Execution(4), costs) do i
        visits[i] += 1
    end
    @test visits == ones(Int, length(visits))
end

@testset "Serial/threaded scan equivalence" begin
    # Create a small PWM model
    weights = Float32[
        0.5 -0.5 0.3
        -0.3 0.7 -0.2
        0.1 0.1 0.8
        -0.2 0.3 -0.1
        -0.3 -0.3 -0.3
    ]
    bg = (0.25f0, 0.25f0, 0.25f0, 0.25f0)
    pwm = PWM("test", weights, bg)

    # Create a batch with mixed-length sequences
    seqs = [
        encode_sequence("ACGTACGTACGTACGT"),
        encode_sequence("ACGTAC"),
        encode_sequence("TTTTGGGGCCCCAAAACGT"),
        encode_sequence("A"),
        encode_sequence("NNNNACGTNNNN"),
    ]
    data = UInt8[]
    offsets = [1]
    for s in seqs
        append!(data, s)
        push!(offsets, length(data) + 1)
    end
    batch = EncodedSequenceBatch(data, offsets)

    # Test all strand policies
    for strands in (ForwardOnly(), ReverseOnly(), BestStrand(), BothStrands())
        serial_result = scan(pwm, batch; strands=strands, execution=Execution())
        for nt in (1, 2, 4)
            threaded_result = scan(pwm, batch; strands=strands, execution=Execution(nt))
            if strands isa BothStrands
                @test threaded_result.forward.data == serial_result.forward.data
                @test threaded_result.forward.offsets == serial_result.forward.offsets
                @test threaded_result.reverse.data == serial_result.reverse.data
                @test threaded_result.reverse.offsets == serial_result.reverse.offsets
            else
                @test threaded_result.data == serial_result.data
                @test threaded_result.offsets == serial_result.offsets
            end
        end
    end
end

@testset "Serial/threaded BaMM scan equivalence" begin
    # Create a small BaMM model (order 1, 3 positions)
    rep = zeros(Float32, 25, 3)  # 5^2 = 25 rows
    for i in 1:25
        rep[i, :] .= Float32(i) / 25
    end
    model = BaMM("test_bamm", rep, 1, 3)

    seqs = [
        encode_sequence("ACGTACGTACGT"),
        encode_sequence("ACG"),
        encode_sequence("TTTTGGGGCCCCAAAACGTAC"),
        encode_sequence("A"),
    ]
    data = UInt8[]
    offsets = [1]
    for s in seqs
        append!(data, s)
        push!(offsets, length(data) + 1)
    end
    batch = EncodedSequenceBatch(data, offsets)

    for strands in (ForwardOnly(), ReverseOnly(), BestStrand(), BothStrands())
        serial_result = scan(model, batch; strands=strands, execution=Execution())
        for nt in (1, 2, 4)
            threaded_result = scan(model, batch; strands=strands, execution=Execution(nt))
            if strands isa BothStrands
                @test threaded_result.forward.data == serial_result.forward.data
                @test threaded_result.reverse.data == serial_result.reverse.data
            else
                @test threaded_result.data == serial_result.data
                @test threaded_result.offsets == serial_result.offsets
            end
        end
    end
end
