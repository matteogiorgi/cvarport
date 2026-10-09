"""
    Constraints(; cap = 1.0)

Fully invested (`sum(w) == 1`), long-only (`w .>= 0`) portfolios with
every weight at most `cap`.
"""
struct Constraints
    cap::Float64
    function Constraints(; cap::Real = 1.0)
        0 < cap <= 1 || throw(ArgumentError("cap must be in (0, 1], got $cap"))
        return new(Float64(cap))
    end
end

function _check_feasible(c::Constraints, N::Integer)
    c.cap * N >= 1 - 1e-12 ||
        throw(ArgumentError("cap = $(c.cap) is infeasible with $N assets"))
    return nothing
end

# Linear programs go to HiGHS (simplex), quadratic ones to Clarabel
# (interior point): HiGHS's active-set QP solver can stall on the tiny
# covariances of daily returns.
function _model(optimizer = HiGHS.Optimizer)
    model = Model(optimizer)
    set_silent(model)
    return model
end

# Sigma scaled to unit average variance. The minimizer does not change,
# but the objective becomes of order one, so the solver's absolute
# tolerances mean the same thing for daily and for annual returns.
_unit(Sigma) = Sigma ./ (tr(Sigma) / size(Sigma, 1))

# Budget, long-only and cap constraints on a fresh weight vector.
function _weights!(model, N::Integer, c::Constraints)
    _check_feasible(c, N)
    @variable(model, 0 <= w[1:N] <= c.cap)
    @constraint(model, sum(w) == 1)
    return w
end

function _solve!(model)
    optimize!(model)
    status = termination_status(model)
    status in (MOI.OPTIMAL, MOI.LOCALLY_SOLVED) || error("solver failed: $status")
    return nothing
end

# Remove solver noise: tiny negative weights, sum slightly off one.
function _clean(w)
    v = max.(w, 0)
    return v ./ sum(v)
end

"""
    min_variance(Sigma, c = Constraints(); mu = nothing, target = nothing) -> w

Minimum-variance portfolio; with `mu` and `target`, the one of minimum
variance among those with expected return `mu' w >= target`.
"""
function min_variance(Sigma::AbstractMatrix, c::Constraints = Constraints();
                      mu = nothing, target = nothing)
    N = size(Sigma, 1)
    model = _model(Clarabel.Optimizer)
    w = _weights!(model, N, c)
    if target !== nothing
        mu === nothing && throw(ArgumentError("a target return needs mu"))
        @constraint(model, dot(mu, w) >= target)
    end
    @objective(model, Min, w' * _unit(Sigma) * w)
    _solve!(model)
    return _clean(value.(w))
end

"""
    max_return(mu, c = Constraints()) -> w

Portfolio of maximum expected return: the right end of both frontiers.
"""
function max_return(mu::AbstractVector, c::Constraints = Constraints())
    model = _model()
    w = _weights!(model, length(mu), c)
    @objective(model, Max, dot(mu, w))
    _solve!(model)
    return _clean(value.(w))
end

"""
    tangency(mu, Sigma, rf, c = Constraints()) -> w

Maximum-Sharpe portfolio, with `rf` in the same units as `mu` (per
period). Solved as the convex problem

    min y' Sigma y   s.t.   (mu - rf)' y = 1,  y >= 0,  y_i <= cap sum_j y_j

and `w = y / sum(y)`, which is equivalent because the Sharpe ratio is
invariant to scaling.
"""
function tangency(mu::AbstractVector, Sigma::AbstractMatrix, rf::Real,
                  c::Constraints = Constraints())
    N = length(mu)
    ex = mu .- rf
    dot(ex, max_return(mu, c)) > 0 ||
        throw(ArgumentError("no feasible portfolio beats the risk-free rate"))
    model = _model(Clarabel.Optimizer)
    @variable(model, y[1:N] >= 0)
    @constraint(model, dot(ex, y) == 1)
    @constraint(model, [i = 1:N], y[i] <= c.cap * sum(y))
    @objective(model, Min, y' * _unit(Sigma) * y)
    _solve!(model)
    return _clean(value.(y))
end

"""
    min_cvar(R, alpha, c = Constraints(); mu = mean of R, target = nothing) -> (w, CVaR)

Minimum historical CVaR at level `alpha` over the scenarios in the rows of
`R` (Rockafellar-Uryasev linear program):

    min  zeta + 1/((1-alpha) S) sum_s u_s
    s.t. u_s >= -r_s' w - zeta,  u_s >= 0,  mu' w >= target

At the optimum zeta is the VaR of the portfolio.
"""
function min_cvar(R::AbstractMatrix, alpha::Real, c::Constraints = Constraints();
                  mu = vec(mean(R, dims = 1)), target = nothing)
    _check_alpha(alpha)
    S, N = size(R)
    model = _model()
    w = _weights!(model, N, c)
    @variable(model, zeta)
    @variable(model, u[1:S] >= 0)
    @constraint(model, [s = 1:S], u[s] >= -sum(R[s, i] * w[i] for i in 1:N) - zeta)
    target === nothing || @constraint(model, dot(mu, w) >= target)
    @objective(model, Min, zeta + sum(u) / ((1 - alpha) * S))
    _solve!(model)
    return _clean(value.(w)), objective_value(model)
end

# Targets from the left end of a frontier to its right end. The last
# point backs off slightly so the solver is not asked for an exact edge.
function _targets(lo, hi, n)
    n >= 2 || throw(ArgumentError("a frontier needs at least two points"))
    hi = max(lo, hi - 1e-9 * max(abs(hi), 1))
    return collect(range(lo, hi; length = n))
end

"""
    mv_frontier(mu, Sigma, c = Constraints(); n = 20) -> (targets, W)

Mean-variance efficient frontier: column `k` of `W` is the
minimum-variance portfolio with expected return at least `targets[k]`.
"""
function mv_frontier(mu::AbstractVector, Sigma::AbstractMatrix, c::Constraints = Constraints();
                     n::Integer = 20)
    lo = dot(mu, min_variance(Sigma, c))
    hi = dot(mu, max_return(mu, c))
    targets = _targets(lo, hi, n)
    W = reduce(hcat, [min_variance(Sigma, c; mu = mu, target = t) for t in targets])
    return targets, W
end

"""
    cvar_frontier(R, alpha, c = Constraints(); mu = mean of R, n = 20) -> (targets, W)

Mean-CVaR efficient frontier, built like [`mv_frontier`](@ref).
"""
function cvar_frontier(R::AbstractMatrix, alpha::Real, c::Constraints = Constraints();
                       mu = vec(mean(R, dims = 1)), n::Integer = 20)
    lo = dot(mu, first(min_cvar(R, alpha, c; mu = mu)))
    hi = dot(mu, max_return(mu, c))
    targets = _targets(lo, hi, n)
    W = reduce(hcat, [first(min_cvar(R, alpha, c; mu = mu, target = t)) for t in targets])
    return targets, W
end
