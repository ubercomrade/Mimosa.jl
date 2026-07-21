using Test
using Mimosa
using SHA

@testset "ExecutionPolicy" begin
    @test SerialExecution() isa ExecutionPolicy
    @test ThreadedExecution() isa ExecutionPolicy
    @test ThreadedExecution(4).ntasks == 4
    @test ThreadedExecution(1).ntasks == 1
    @test_throws ArgumentError ThreadedExecution(0)

    # show
    @test occursin("SerialExecution", sprint(show, SerialExecution()))
    @test occursin("ThreadedExecution", sprint(show, ThreadedExecution(2)))
end

@testset "CLI execution policy" begin
    parsed = Mimosa.CLIParsed("profile")
    parsed.options["threads"] = string(Threads.nthreads())
    policy = Mimosa._execution_policy(parsed)
    if Threads.nthreads() == 1
        @test policy isa SerialExecution
    else
        @test policy == ThreadedExecution(Threads.nthreads())
    end

    parsed.options["threads"] = string(Threads.nthreads() + 1)
    @test_throws Mimosa.CLIError Mimosa._execution_policy(parsed)
end

@testset "_parallel_for serial" begin
    # Basic serial execution
    results = Vector{Int}(undef, 10)
    Mimosa._parallel_for(SerialExecution(), 10) do i
        results[i] = i * 2
    end
    @test results == collect(2:2:20)
end

@testset "_parallel_for threaded" begin
    # Threaded execution should produce same results as serial
    for ntasks in (1, 2, 4)
        results_t = Vector{Int}(undef, 100)
        Mimosa._parallel_for(ThreadedExecution(ntasks), 100) do i
            results_t[i] = i^2
        end
        results_s = Vector{Int}(undef, 100)
        Mimosa._parallel_for(SerialExecution(), 100) do i
            results_s[i] = i^2
        end
        @test results_t == results_s
    end

    # Edge case: n=1
    results = Vector{Int}(undef, 1)
    Mimosa._parallel_for(ThreadedExecution(4), 1) do i
        results[i] = 42
    end
    @test results == [42]

    # Edge case: n=0
    Mimosa._parallel_for(ThreadedExecution(4), 0) do i
        @test false  # should not execute
    end
end

@testset "bounded and weighted execution" begin
    @test Mimosa._effective_ntasks(ThreadedExecution(100), 100) == Threads.nthreads()

    visits = zeros(Int, 37)
    costs = [i % 7 == 0 ? 1000 : 1 for i in eachindex(visits)]
    Mimosa._parallel_for_weighted(ThreadedExecution(4), costs) do i
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
        serial_result = scan(pwm, batch; strands=strands, execution=SerialExecution())
        for nt in (1, 2, 4)
            threaded_result = scan(
                pwm, batch; strands=strands, execution=ThreadedExecution(nt)
            )
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
        serial_result = scan(model, batch; strands=strands, execution=SerialExecution())
        for nt in (1, 2, 4)
            threaded_result = scan(
                model, batch; strands=strands, execution=ThreadedExecution(nt)
            )
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
