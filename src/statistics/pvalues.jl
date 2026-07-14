# P-value, E-value, and Benjamini-Hochberg FDR adjustment.

"""
    BenjaminiHochberg

Type marker for the Benjamini-Hochberg FDR adjustment method.
"""
struct BenjaminiHochberg end

"""
    adjusted_pvalues(pvalues; method=BenjaminiHochberg())

Compute FDR-adjusted p-values using the Benjamini-Hochberg procedure.

The BH procedure:
1. Sort p-values: p_(1) ≤ p_(2) ≤ ... ≤ p_(m)
2. Starting from the largest, compute q_(i) = min(q_(i+1), p_(i) * m / i)
3. Cap at 1.0
4. Return in original order
"""
function adjusted_pvalues(
    pvalues::AbstractVector{<:Real}; method::BenjaminiHochberg=BenjaminiHochberg()
)
    n = length(pvalues)
    n == 0 && return Float64[]

    p = Float64.(collect(pvalues))
    all(isfinite, p) && all((0.0 .<= p) .& (p .<= 1.0)) ||
        throw(ArgumentError("p-values must be finite and lie in [0, 1]."))
    order = sortperm(p)
    sorted_p = p[order]

    adj = Vector{Float64}(undef, n)
    # Start from the largest p-value
    adj[n] = min(sorted_p[n] * n / n, 1.0)

    @inbounds for i in (n - 1):-1:1
        val = sorted_p[i] * n / i
        adj[i] = min(adj[i + 1], val)
    end

    # Cap at 1.0
    @inbounds for i in 1:n
        if adj[i] > 1.0
            adj[i] = 1.0
        end
    end

    # Restore original order
    result = Vector{Float64}(undef, n)
    @inbounds for j in 1:n
        result[order[j]] = adj[j]
    end

    return result
end

"""
    evalue(pvalue::Real, effective_n::Int)

Compute E-value: `p * effective_n`. The effective number of targets accounts
for multiple testing in the comparison set.
"""
function evalue(pvalue::Real, effective_n::Int)
    p = Float64(pvalue)
    isfinite(p) && 0.0 <= p <= 1.0 ||
        throw(ArgumentError("p-value must be finite and lie in [0, 1]."))
    effective_n >= 0 || throw(ArgumentError("effective_n must be non-negative."))
    return p * effective_n
end

function evalue(pvalue::Real, effective_n::Real)
    isfinite(effective_n) && effective_n >= 0 ||
        throw(ArgumentError("effective_n must be finite and non-negative."))
    isinteger(effective_n) || throw(ArgumentError("effective_n must be an integer."))
    return evalue(pvalue, Int(effective_n))
end

"""
    pvalue(gev::GEVFit, score::Real)

Upper-tail p-value for a score under the fitted GEV distribution.
Equivalent to [`survival`](@ref).
"""
function pvalue(gev::GEVFit, score::Real)
    return survival(gev, score)
end
