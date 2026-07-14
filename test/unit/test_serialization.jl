using Test
using Mimosa

@testset "to_dict produces expected keys" begin
    result = ComparisonResult("q", "t", 1.0f0, 0, "++", "co")
    d = Mimosa.to_dict(result)
    @test Set(keys(d)) ==
        Set(["query", "target", "score", "offset", "orientation", "metric"])
    @test d["score"] == 1.0
    @test d["offset"] == 0
end

@testset "to_json round-trips query/target" begin
    result = ComparisonResult("pwm_model", "pwm_model", 1.0f0, 0, "++", "co")
    s = Mimosa.to_json(result)
    @test contains(s, "\"query\"")
    @test contains(s, "\"pwm_model\"")
    @test contains(s, "\"orientation\"")
    @test contains(s, "\"++\"")
    @test contains(s, "\"co\"")
end

@testset "annotated result declares its JSON contract version" begin
    result = ComparisonResult("q", "t", 1.0f0, 0, "++", "co")
    annotated = AnnotatedResult(result; p_value=0.1, adj_p_value=0.1, e_value=0.1)
    payload = Mimosa.to_dict(annotated)
    @test payload["annotation_schema_version"] == ANNOTATED_RESULT_SCHEMA_VERSION
end
