# Native GEV fitting and survival function.
#
# Uses the textbook shape-parameter convention `k` (sign flipped from SciPy's `c`).
# See ADR 0005 for the parameterization audit and compatibility plan.
#
# GEV CDF (textbook convention, shape k, location μ, scale σ > 0):
#   k ≠ 0: F(x) = exp(-(1 + k*(x-μ)/σ)^(-1/k))   where 1 + k*(x-μ)/σ > 0
#   k = 0: F(x) = exp(-exp(-(x-μ)/σ))
#
# Support:
#   k > 0 (Fréchet): (μ - σ/k, +∞)
#   k = 0 (Gumbel):  (-∞, +∞)
#   k < 0 (Weibull): (-∞, μ - σ/k)

"""
    GEVFit

Result of a successful GEV MLE fit.

Fields:
- `shape::Float64`: shape parameter `k` (textbook convention, sign flipped from SciPy's `c`).
- `location::Float64`: location parameter `μ`.
- `scale::Float64`: scale parameter `σ` (always positive).
- `converged::Bool`: whether the optimizer converged.
- `iterations::Int`: number of iterations used.
- `loglikelihood::Float64`: log-likelihood at the fitted parameters.
"""
struct GEVFit
    shape::Float64
    location::Float64
    scale::Float64
    converged::Bool
    iterations::Int
    loglikelihood::Float64

    function GEVFit(
        shape::Float64,
        location::Float64,
        scale::Float64,
        converged::Bool,
        iterations::Int,
        loglikelihood::Float64,
    )
        all(isfinite, (shape, location, scale, loglikelihood)) ||
            throw(ArgumentError("GEV parameters must be finite."))
        scale > 0 || throw(ArgumentError("GEV scale must be positive."))
        iterations >= 0 || throw(ArgumentError("GEV iterations must be non-negative."))
        return new(shape, location, scale, converged, iterations, loglikelihood)
    end
end

function GEVFit(
    shape::Real,
    location::Real,
    scale::Real,
    converged::Bool,
    iterations::Int,
    loglikelihood::Real,
)
    return GEVFit(
        Float64(shape),
        Float64(location),
        Float64(scale),
        converged,
        iterations,
        Float64(loglikelihood),
    )
end

"""
    GEVFitFailure

Typed failure result when GEV fitting cannot produce a valid distribution.

Fields:
- `message::String`: human-readable explanation.
- `n_scores::Int`: number of input scores.
- `iterations::Int`: iterations spent before failure.
"""
struct GEVFitFailure
    message::String
    n_scores::Int
    iterations::Int
end

"""
    GEVFitResult

Union of [`GEVFit`](@ref) and [`GEVFitFailure`](@ref).
"""
const GEVFitResult = Union{GEVFit,GEVFitFailure}

# ---------------------------------------------------------------------------
# Negative log-likelihood
# ---------------------------------------------------------------------------

"""
    gev_nll(params, data)

Negative log-likelihood of GEV distribution at `params = [k, μ, log_σ]` given
`data`. Returns `Inf` when the support constraint `1 + k*(x-μ)/σ > 0` is
violated for any data point.
"""
function gev_nll(params::AbstractVector{Float64}, data::Vector{Float64})
    k = params[1]
    mu = params[2]
    log_sigma = params[3]
    sigma = exp(log_sigma)
    n = length(data)

    if abs(k) < 1e-12
        # Gumbel limit: NLL = n*log(σ) + Σ[z_i + exp(-z_i)]
        nll = n * log_sigma
        @inbounds for i in 1:n
            z = (data[i] - mu) / sigma
            nll += z + exp(-z)
        end
        return nll
    end

    inv_k = 1.0 / k
    one_plus_inv_k = 1.0 + inv_k
    nll = n * log_sigma

    @inbounds for i in 1:n
        z = (data[i] - mu) / sigma
        y = 1.0 + k * z
        if y <= 0.0
            return Inf
        end
        nll += one_plus_inv_k * log(y) + y^(-inv_k)
    end
    return nll
end

# ---------------------------------------------------------------------------
# Numerical gradient (central differences with fallback)
# ---------------------------------------------------------------------------

_vec_dot(a::Vector{Float64}, b::Vector{Float64}) = sum(a[i] * b[i] for i in eachindex(a))

function _numerical_gradient(f, x::Vector{Float64}; h::Float64=1e-5)
    n = length(x)
    g = Vector{Float64}(undef, n)
    # Pre-allocate work vectors to avoid per-dimension allocations
    xp = copy(x)
    xm = copy(x)
    for i in 1:n
        xp[i] = x[i] + h
        fp = f(xp)
        xm[i] = x[i] - h
        fm = f(xm)
        # Restore original value
        xp[i] = x[i]
        xm[i] = x[i]
        if isinf(fp) && isinf(fm)
            g[i] = 0.0
        elseif isinf(fp)
            f0 = f(x)
            g[i] = (f0 - fm) / h
        elseif isinf(fm)
            f0 = f(x)
            g[i] = (fp - f0) / h
        else
            g[i] = (fp - fm) / (2.0 * h)
        end
    end
    return g
end

function _numerical_gradient!(g::Vector{Float64}, f, x::Vector{Float64}; h::Float64=1e-5)
    length(g) == length(x) || throw(ArgumentError("gradient buffer has wrong length."))
    n = length(x)
    xp = copy(x)
    xm = copy(x)
    for i in 1:n
        xp[i] = x[i] + h
        fp = f(xp)
        xm[i] = x[i] - h
        fm = f(xm)
        xp[i] = x[i]
        xm[i] = x[i]
        if isinf(fp) && isinf(fm)
            g[i] = 0.0
        elseif isinf(fp)
            g[i] = (f(x) - fm) / h
        elseif isinf(fm)
            g[i] = (fp - f(x)) / h
        else
            g[i] = (fp - fm) / (2.0 * h)
        end
    end
    return g
end

# ---------------------------------------------------------------------------
# BFGS optimizer with backtracking line search
# ---------------------------------------------------------------------------

struct BFGSResult
    x::Vector{Float64}
    converged::Bool
    iterations::Int
    gradient_norm::Float64
end

function _make_identity(n::Int)
    H = Matrix{Float64}(undef, n, n)
    for i in 1:n, j in 1:n
        H[i, j] = i == j ? 1.0 : 0.0
    end
    return H
end

function _bfgs_optimize(f, x0::Vector{Float64}; max_iter::Int=500, tol::Float64=1e-8)
    n = length(x0)
    x = copy(x0)
    H = _make_identity(n)
    g = _numerical_gradient(f, x)

    # Pre-allocate reusable work vectors
    p = Vector{Float64}(undef, n)
    x_new = Vector{Float64}(undef, n)
    g_new = Vector{Float64}(undef, n)
    s = Vector{Float64}(undef, n)
    y_vec = Vector{Float64}(undef, n)
    Hy = Vector{Float64}(undef, n)

    for iter in 1:max_iter
        gnorm = sqrt(sum(abs2, g))
        if gnorm < tol
            return BFGSResult(x, true, iter, gnorm)
        end

        # Search direction: p = -H * g
        for i in 1:n
            p[i] = -sum(H[i, j] * g[j] for j in 1:n)
        end
        pnorm = sqrt(sum(abs2, p))
        if pnorm < 1e-15
            return BFGSResult(x, true, iter, gnorm)
        end

        # Backtracking line search with Armijo condition
        f0 = f(x)
        dg = _vec_dot(g, p)
        alpha = 1.0
        for i in 1:n
            x_new[i] = x[i] + alpha * p[i]
        end
        f_new = f(x_new)

        backtrack = 0
        while f_new > f0 + 1e-4 * alpha * dg && backtrack < 50
            alpha *= 0.5
            for i in 1:n
                x_new[i] = x[i] + alpha * p[i]
            end
            f_new = f(x_new)
            backtrack += 1
        end

        if backtrack >= 50 || isinf(f_new)
            # Line search failed — try a small step
            alpha = 1e-8
            for i in 1:n
                x_new[i] = x[i] + alpha * p[i]
            end
            f_new = f(x_new)
            if isinf(f_new) || f_new >= f0
                return BFGSResult(x, false, iter, gnorm)
            end
        end

        _numerical_gradient!(g_new, f, x_new)
        for i in 1:n
            s[i] = x_new[i] - x[i]
            y_vec[i] = g_new[i] - g[i]
        end
        sy = _vec_dot(s, y_vec)

        if sy > 1e-12
            rho = 1.0 / sy
            # BFGS update: H = (I - ρ*s*y') * H * (I - ρ*y*s') + ρ*s*s'
            for i in 1:n
                Hy[i] = sum(H[i, j] * y_vec[j] for j in 1:n)
            end
            yHy = _vec_dot(y_vec, Hy)
            coef = rho + rho * rho * yHy
            for i in 1:n, j in 1:n
                H[i, j] += coef * s[i] * s[j] - rho * (Hy[i] * s[j] + s[i] * Hy[j])
            end
        end

        for i in 1:n
            x[i] = x_new[i]
            g[i] = g_new[i]
        end
    end

    gnorm = sqrt(sum(abs2, g))
    return BFGSResult(x, gnorm < tol * 100, max_iter, gnorm)
end

# ---------------------------------------------------------------------------
# Public fit API
# ---------------------------------------------------------------------------

"""
    fit_gev(scores; max_iter=500, tol=1e-8) -> GEVFitResult

Fit a Generalized Extreme Value distribution to `scores` via maximum
likelihood. Returns a [`GEVFit`](@ref) on success or a [`GEVFitFailure`](@ref)
on failure (degenerate sample, non-convergence, non-finite parameters).

The shape parameter uses the textbook convention `k` (sign flipped from
SciPy's `c`): `k = -c`.
"""
function fit_gev(scores::AbstractVector{<:Real}; max_iter::Int=500, tol::Float64=1e-8)
    max_iter > 0 || throw(ArgumentError("max_iter must be positive."))
    isfinite(tol) && tol > 0 || throw(ArgumentError("tol must be finite and positive."))
    data = sort!(Float64.(collect(scores)))
    n = length(data)

    if n < 3
        return GEVFitFailure("GEV fit requires at least 3 scores (got $n).", n, 0)
    end

    # Check for degenerate (constant) sample
    if data[end] == data[1]
        return GEVFitFailure(
            "Degenerate sample: all $n scores are identical ($(data[1])).", n, 0
        )
    end

    # Check for NaN/Inf
    if any(!isfinite, data)
        return GEVFitFailure("Sample contains non-finite values.", n, 0)
    end

    # Method-of-moments initialization (Gumbel estimates)
    mu0 = sum(data) / n
    var0 = sum((data .- mu0) .^ 2) / (n - 1)
    sigma0 = sqrt(max(var0, 1e-10)) * sqrt(6.0) / pi
    if sigma0 < 1e-10
        sigma0 = 1e-6
    end

    # Parameterize: [k, μ, log_σ]  (log_σ ensures σ > 0)
    x0 = [0.0, mu0, log(sigma0)]

    nll = x -> gev_nll(x, data)
    result = _bfgs_optimize(nll, x0; max_iter=max_iter, tol=tol)

    k = result.x[1]
    mu = result.x[2]
    sigma = exp(result.x[3])

    # Validate fitted parameters
    if !isfinite(k) || !isfinite(mu) || !isfinite(sigma) || sigma <= 0
        return GEVFitFailure(
            "Fitted parameters are non-finite or scale ≤ 0: k=$k, μ=$mu, σ=$sigma.",
            n,
            result.iterations,
        )
    end

    # Check support constraint at all data points
    if abs(k) > 1e-12
        for x in data
            y = 1.0 + k * (x - mu) / sigma
            if y <= 0
                return GEVFitFailure(
                    "Support constraint violated at data points.", n, result.iterations
                )
            end
        end
    end

    ll = -nll(result.x)

    return GEVFit(k, mu, sigma, result.converged, result.iterations, ll)
end

# ---------------------------------------------------------------------------
# Survival function (upper-tail p-value)
# ---------------------------------------------------------------------------

"""
    survival(gev::GEVFit, x::Real)

Upper-tail survival probability `P(X > x)` for a fitted GEV distribution.
Clamped to `[0, 1]`.

For the Gumbel case (|k| < 1e-12), uses `-expm1(-exp(-z))` for precision.
For general k, uses `-expm1(-y^(-1/k))` where `y = 1 + k*z`.
"""
function survival(gev::GEVFit, x::Real)
    k = gev.shape
    mu = gev.location
    sigma = gev.scale
    xv = Float64(x)

    if abs(k) < 1e-12
        # Gumbel: SF = 1 - exp(-exp(-z)) = -expm1(-exp(-z))
        z = (xv - mu) / sigma
        return clamp(-expm1(-exp(-z)), 0.0, 1.0)
    end

    z = (xv - mu) / sigma
    y = 1.0 + k * z

    if y <= 0.0
        # Outside support
        if k > 0
            # Fréchet: lower-bounded, x below lower bound → SF = 1
            return 1.0
        else
            # Weibull: upper-bounded, x above upper bound → SF = 0
            return 0.0
        end
    end

    t = y^(-1.0 / k)
    return clamp(-expm1(-t), 0.0, 1.0)
end

"""
    cdf(gev::GEVFit, x::Real)

CDF `P(X ≤ x)` for a fitted GEV distribution.
"""
function cdf(gev::GEVFit, x::Real)
    return 1.0 - survival(gev, x)
end

"""
    scipy_params(gev::GEVFit)

Return `(c, loc, scale)` in SciPy's convention (c = -k).
Useful for comparison with frozen oracle fixtures.
"""
scipy_params(gev::GEVFit) = (-gev.shape, gev.location, gev.scale)
