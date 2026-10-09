# The synthetic example through the library API instead of the command
# line: same analysis, but the results go to CSV files that the R
# scripts in scripts/ turn into the figures of the README.
#
#   julia --project=. examples/synthetic.jl > examples/prices.csv
#   julia --project=. examples/analysis.jl
#
# Writes examples/portfolios.csv, examples/contributions.csv and
# examples/frontier.csv.

using CVaRPort, Printf

const DIR = @__DIR__
const ALPHA, RF, CAP = 0.95, 0.02, 0.7

_, tickers, P = read_prices(joinpath(DIR, "prices.csv"))
w = read_weights(joinpath(DIR, "weights.txt"), tickers)
R = simple_returns(P)
res = analyze(R, w; alpha = ALPHA, rf = RF, cap = CAP)

open(joinpath(DIR, "portfolios.csv"), "w") do io
    println(io, join(["portfolio", "mean", "vol", "sharpe", "var", "cvar", tickers...], ','))
    for p in res.portfolios
        s = p.stats
        name = "\"" * p.name * "\""                   # some names contain commas
        println(io, join([name; string.([s.mean, s.vol, s.sharpe, s.var, s.cvar]); string.(p.w)], ','))
    end
end

open(joinpath(DIR, "contributions.csv"), "w") do io
    vol, cvar = res.contributions.vol, res.contributions.cvar
    println(io, "asset,weight,vol,cvar")
    for (i, t) in enumerate(tickers)
        @printf(io, "%s,%.6f,%.6f,%.6f\n", t, w[i], vol[i] / sum(vol), cvar[i] / sum(cvar))
    end
end

write_frontier(joinpath(DIR, "frontier.csv"), tickers, R; alpha = ALPHA, cap = CAP)
