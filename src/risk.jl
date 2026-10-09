# Losses are L = -R w: positive numbers are money lost per unit invested.

# Index of the alpha-quantile in the sorted losses. The tolerance protects
# against alpha * S landing just above an integer in floating point.
_tail_index(S::Integer, alpha::Real) = clamp(ceil(Int, alpha * S - 1e-9), 1, S)

_check_alpha(alpha) = 0 < alpha < 1 || throw(ArgumentError("alpha must be in (0, 1), got $alpha"))

"""
    var_cvar(L, alpha) -> (VaR, CVaR)
    var_cvar(R, w, alpha) -> (VaR, CVaR)

Historical VaR and CVaR of the losses `L` (or of `-R w`), with the
discrete Rockafellar-Uryasev definition

    CVaR = min_zeta  zeta + 1/((1-alpha) S) sum_s max(L_s - zeta, 0)

whose minimizer zeta is the VaR. This is exactly the quantity minimized
by [`min_cvar`](@ref), so the two always agree.
"""
function var_cvar(L::AbstractVector{<:Real}, alpha::Real)
    _check_alpha(alpha)
    S = length(L)
    S > 0 || throw(ArgumentError("no scenarios"))
    zeta = partialsort(L, _tail_index(S, alpha))
    cvar = zeta + sum(l -> max(l - zeta, 0), L) / ((1 - alpha) * S)
    return zeta, cvar
end

var_cvar(R::AbstractMatrix, w::AbstractVector, alpha::Real) = var_cvar(-(R * w), alpha)

"""
    cvar_contributions(R, w, alpha) -> c

Euler decomposition of the historical CVaR: `c[i]` is the part of the
CVaR due to asset `i`, and `sum(c) == CVaR` exactly.
"""
function cvar_contributions(R::AbstractMatrix, w::AbstractVector, alpha::Real)
    _check_alpha(alpha)
    S = size(R, 1)
    L = -(R * w)
    k = _tail_index(S, alpha)
    sstar = sortperm(L)[k]
    zeta = L[sstar]
    tail = L .> zeta
    m = (1 - alpha) * S
    g = (vec(sum(-R[tail, :], dims = 1)) .- R[sstar, :] .* (m - count(tail))) ./ m
    return w .* g
end

"""
    vol_contributions(Sigma, w) -> c

Euler decomposition of the volatility: `sum(c) == sqrt(w' Sigma w)`.
"""
function vol_contributions(Sigma::AbstractMatrix, w::AbstractVector)
    sigma = sqrt(dot(w, Sigma, w))
    return w .* (Sigma * w) ./ sigma
end

"""
    portfolio_stats(mu, Sigma, R, w; alpha, rf, periods) -> NamedTuple

Mean and volatility annualized with `periods` per year (linear
scaling of simple returns), Sharpe ratio against the annual risk-free
rate `rf`, and one-period historical VaR and CVaR.
"""
function portfolio_stats(mu, Sigma, R, w; alpha = 0.95, rf = 0.0, periods = 252)
    mean_a = periods * dot(mu, w)
    vol_a = sqrt(periods * dot(w, Sigma, w))
    var, cvar = var_cvar(R, w, alpha)
    return (mean = mean_a, vol = vol_a, sharpe = (mean_a - rf) / vol_a, var = var, cvar = cvar)
end
