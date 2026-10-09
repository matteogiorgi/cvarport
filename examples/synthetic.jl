# Synthetic daily prices to try cvarport without a network connection:
#
#   julia --project=. examples/synthetic.jl > examples/prices.csv
#   julia --project=. bin/cvarport.jl examples/prices.csv examples/weights.txt --rf 0.02 --cap 0.5
#
# Four fake assets (equity-like, bond-like, gold-like, small-cap-like),
# about five years of trading days, fat-tailed correlated returns.

using LinearAlgebra, Random

rng = MersenneTwister(2026)
tickers = ["EQTY", "BOND", "GOLD", "SMCP"]
mu = [0.07, 0.02, 0.04, 0.09] ./ 252                # annual drift / 252
vol = [0.16, 0.05, 0.14, 0.24] ./ sqrt(252)
rho = [1.0 0.1 0.0 0.8
       0.1 1.0 0.2 0.0
       0.0 0.2 1.0 0.0
       0.8 0.0 0.0 1.0]
C = cholesky(Symmetric(vol .* rho .* vol')).L

T = 1260
P = fill(100.0, 1, 4)
for _ in 1:T
    z = randn(rng, 4) ./ sqrt(sum(abs2, randn(rng, 4)) / 4) ./ sqrt(2)   # Student t, 4 dof
    global P = vcat(P, P[end:end, :] .* (1 .+ (mu .+ C * z))')
end

println(join(["date"; tickers], ','))
for (t, row) in enumerate(eachrow(P))
    println(join([string("d", lpad(t, 5, '0')); string.(round.(row; digits = 4))], ','))
end
