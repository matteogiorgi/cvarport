"""
    analyze(R, w; alpha = 0.95, rf = 0.0, periods = 252, cap = 1.0, shrink = true) -> NamedTuple

Compare the portfolio `w` with the reference portfolios computed on
the same returns `R` and constraints. `rf` is the annual risk-free
rate; `periods` the number of return periods in a year. With `shrink`
the covariance is Ledoit-Wolf, otherwise the sample one.

The comparison portfolios are: minimum variance, tangency (when some
portfolio beats `rf`), minimum CVaR, and the minimum-variance and
minimum-CVaR portfolios with expected return at least that of `w`.
"""
function analyze(R::AbstractMatrix, w::AbstractVector; alpha::Real = 0.95, rf::Real = 0.0,
                 periods::Integer = 252, cap::Real = 1.0, shrink::Bool = true)
    size(R, 2) == length(w) || throw(DimensionMismatch("R has $(size(R, 2)) assets, w has $(length(w))"))
    c = Constraints(cap = cap)
    mu = vec(mean(R, dims = 1))
    Sigma, delta = shrink ? ledoit_wolf(R) : (cov(R, corrected = false), 0.0)
    within_cap = all(<=(cap + 1e-9), w)
    # If w breaks the cap its return may be out of reach: clamp the target.
    target = min(dot(mu, w), dot(mu, max_return(mu, c)))

    named = Pair{String,Vector{Float64}}["current" => w, "min variance" => min_variance(Sigma, c)]
    try
        push!(named, "tangency" => tangency(mu, Sigma, rf / periods, c))
    catch err
        err isa ArgumentError || rethrow()
    end
    push!(named, "min CVaR" => first(min_cvar(R, alpha, c; mu = mu)))
    push!(named, "min var, mean >= current" => min_variance(Sigma, c; mu = mu, target = target))
    push!(named, "min CVaR, mean >= current" => first(min_cvar(R, alpha, c; mu = mu, target = target)))

    portfolios = [(name = name, w = v,
                   stats = portfolio_stats(mu, Sigma, R, v; alpha = alpha, rf = rf, periods = periods))
                  for (name, v) in named]
    contributions = (vol = vol_contributions(Sigma, w) .* sqrt(periods),
                     cvar = cvar_contributions(R, w, alpha))
    return (; mu, Sigma, delta, alpha, rf, periods, cap, within_cap, scenarios = size(R, 1),
            portfolios, contributions)
end

_pct(x) = @sprintf("%7.2f%%", 100x)

"""
    print_report([io], tickers, result)

Human-readable report of an [`analyze`](@ref) result.
"""
function print_report(io::IO, tickers::AbstractVector{<:AbstractString}, res)
    alpha = res.alpha
    lvl = round(Int, 100 * alpha)
    @printf(io, "cvarport: %d return periods, %d assets, alpha = %.3g, cap = %.2f, Ledoit-Wolf delta = %.3f\n",
            res.scenarios, length(tickers), alpha, res.cap, res.delta)
    @printf(io, "mean and vol annualized (x%d), rf = %.2f%%; VaR and CVaR are one-period losses\n\n",
            res.periods, 100res.rf)
    res.within_cap || println(io, "warning: current weights exceed the cap\n")

    names = [first(t, 8) for t in tickers]
    @printf(io, "%-26s %8s %8s %7s %8s %8s  |", "portfolio", "mean", "vol", "Sharpe",
            "VaR$lvl", "CVaR$lvl")
    foreach(t -> @printf(io, " %8s", t), names)
    println(io)
    for p in res.portfolios
        s = p.stats
        @printf(io, "%-26s %s %s %7.3f %s %s  |", p.name, _pct(s.mean), _pct(s.vol),
                s.sharpe, _pct(s.var), _pct(s.cvar))
        foreach(x -> @printf(io, " %8.3f", x), p.w)
        println(io)
    end

    cur = res.portfolios[1]
    vol, cv = res.contributions.vol, res.contributions.cvar
    println(io, "\nrisk contributions of the current portfolio (shares of the total)")
    @printf(io, "%-8s %8s %8s %8s\n", "asset", "weight", "vol", "CVaR$lvl")
    for (i, t) in enumerate(names)
        cur.w[i] > 0 || continue
        @printf(io, "%-8s %8.3f %s %s\n", t, cur.w[i], _pct(vol[i] / sum(vol)), _pct(cv[i] / sum(cv)))
    end

    best = res.portfolios[end]
    dcvar = best.stats.cvar - cur.stats.cvar
    @printf(io, "\nwith expected return at least the current one, CVaR%d goes from %.2f%% to %.2f%% (%+.1f%% relative)\n",
            lvl, 100cur.stats.cvar, 100best.stats.cvar, 100dcvar / cur.stats.cvar)
    return nothing
end

print_report(tickers, res) = print_report(stdout, tickers, res)

"""
    write_frontier(path, tickers, R; alpha = 0.95, periods = 252, cap = 1.0, shrink = true, n = 20)

Write both efficient frontiers to a CSV with columns
`frontier,mean,vol,cvar,<tickers>...` (mean and vol annualized), for
plotting with any tool.
"""
function write_frontier(path::AbstractString, tickers, R::AbstractMatrix; alpha::Real = 0.95,
                        periods::Integer = 252, cap::Real = 1.0, shrink::Bool = true,
                        n::Integer = 20)
    c = Constraints(cap = cap)
    mu = vec(mean(R, dims = 1))
    Sigma = shrink ? first(ledoit_wolf(R)) : cov(R, corrected = false)
    open(path, "w") do io
        println(io, join(["frontier", "mean", "vol", "cvar", tickers...], ','))
        for (label, (_, W)) in ("mv" => mv_frontier(mu, Sigma, c; n = n),
                                "cvar" => cvar_frontier(R, alpha, c; mu = mu, n = n))
            for w in eachcol(W)
                s = portfolio_stats(mu, Sigma, R, w; alpha = alpha, periods = periods)
                println(io, join([label; string.([s.mean, s.vol, s.cvar]); string.(w)], ','))
            end
        end
    end
    return path
end
