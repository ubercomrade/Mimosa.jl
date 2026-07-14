using Test
using Mimosa
using JET

@testset "JET prepared-profile entry point" begin
    profile = ScoreProfile("profile", build_ragged([Float32[0.1, 0.8, 0.2, 0.7]]))
    target = ScoreProfile("target", build_ragged([Float32[0.2, 0.7, 0.3, 0.6]]))
    prepared_profile = prepare_profile(profile)
    prepared_target = prepare_profile(target)
    @test_opt compare(prepared_profile, prepared_target; metric=OverlapCoefficient())
end
