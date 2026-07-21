using Test
using Mimosa

@testset "profile type stability smoke" begin
    p1 = ScoreProfile("q", build_ragged([Float32[0.1, 0.8, 0.2, 0.7]]))
    p2 = ScoreProfile("t", build_ragged([Float32[0.2, 0.7, 0.3, 0.6]]))
    result = compare(p1, p2; metric=:co)
    @test result isa ComparisonResult
    @test result.metric == "co"

    scores = Float32[sin(index) for index in 1:1_000]
    table = @inferred fit(
        HybridEmpiricalLogTail(256), scores; tail_logfpr=2.0, execution=ThreadedExecution(2)
    )
    @test table isa HybridLogTailTable
end
