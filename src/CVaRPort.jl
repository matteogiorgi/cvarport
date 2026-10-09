"""
    CVaRPort

Risk calculator for a given long-only portfolio: historical VaR/CVaR,
risk contributions, and comparison against the mean-variance and
mean-CVaR efficient frontiers.
"""
module CVaRPort

using DelimitedFiles, LinearAlgebra, Printf, Statistics
using Clarabel, HiGHS, JuMP

export read_prices, read_weights, simple_returns
export ledoit_wolf
export var_cvar, cvar_contributions, vol_contributions, portfolio_stats
export Constraints, min_variance, max_return, tangency, min_cvar
export mv_frontier, cvar_frontier
export analyze, print_report, write_frontier

include("data.jl")
include("estimate.jl")
include("risk.jl")
include("optim.jl")
include("report.jl")

end
