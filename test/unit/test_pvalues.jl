using Test
using Mimosa

@testset "BH FDR adjustment" begin
    # Standard example: 5 p-values
    pvals = [0.01, 0.02, 0.03, 0.04, 0.05]
    adj = adjusted_pvalues(pvals)

    @test length(adj) == 5
    # All adjusted p-values should be >= original p-values
    @test all(adj .>= pvals .- 1e-10)
    # All should be <= 1.0
    @test all(adj .<= 1.0)
    # Monotonically non-decreasing in sorted order
    sorted_adj = adj[sortperm(pvals)]
    for i in 1:(length(sorted_adj) - 1)
        @test sorted_adj[i] <= sorted_adj[i + 1] + 1e-10
    end

    # Single p-value: adjusted = p
    @test adjusted_pvalues([0.05]) == [0.05]

    # Empty: returns empty
    @test isempty(adjusted_pvalues(Float64[]))
    @test_throws ArgumentError adjusted_pvalues([NaN])
    @test_throws ArgumentError adjusted_pvalues([-0.1])
    @test_throws ArgumentError adjusted_pvalues([1.1])

    # All p-values = 0: adjusted = 0
    @test all(adjusted_pvalues([0.0, 0.0, 0.0]) .== 0.0)

    # All p-values = 1: adjusted = 1
    @test all(adjusted_pvalues([1.0, 1.0, 1.0]) .== 1.0)

    # Order preservation
    pvals = [0.05, 0.01, 0.04, 0.02, 0.03]
    adj = adjusted_pvalues(pvals)
    order = sortperm(pvals)
    sorted_adj = adj[order]
    for i in 1:(length(sorted_adj) - 1)
        @test sorted_adj[i] <= sorted_adj[i + 1] + 1e-10
    end
end

@testset "E-value" begin
    @test evalue(0.05, 100) ≈ 5.0
    @test evalue(0.0, 100) ≈ 0.0
    @test evalue(1.0, 50) ≈ 50.0
    @test evalue(0.001, 1000) ≈ 1.0
    @test_throws ArgumentError evalue(-0.1, 10)
    @test_throws ArgumentError evalue(0.1, -1)
end

@testset "pvalue from GEV" begin
    gev = GEVFit(0.0, 0.0, 1.0, true, 1, 0.0)
    pv = pvalue(gev, 0.0)
    @test pv ≈ 1.0 - exp(-exp(0.0)) atol = 1e-10
    @test 0.0 ≤ pv ≤ 1.0
end
