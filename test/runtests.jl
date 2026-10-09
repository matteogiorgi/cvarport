using CVaRPort, LinearAlgebra, Random, Statistics, Test

# Scenarios with correlated, fat-tailed returns (Student t, 4 dof).
function scenarios(S, N; seed = 1)
    rng = MersenneTwister(seed)
    A = 0.006 .* randn(rng, N, N)
    C = cholesky(Symmetric(A * A' + 1e-5I)).L
    mu = range(0.0002, 0.0006; length = N)
    chi2 = [sum(abs2, randn(rng, 4)) for _ in 1:S]
    Z = randn(rng, S, N) ./ sqrt.(chi2 ./ 4) ./ sqrt(2)
    return mu' .+ Z * C'
end

const Sigma3 = [0.040 0.006 0.004
                0.006 0.090 0.010
                0.004 0.010 0.0625]

@testset "CVaRPort" begin
    @testset "var_cvar by hand" begin
        L = collect(1.0:10.0)
        @test all(isapprox.(var_cvar(L, 0.8), (8.0, 9.5)))  # mean of the worst 2
        zeta, cv = var_cvar(L, 0.75)                     # worst 2.5 scenarios
        @test zeta == 8.0
        @test isapprox(cv, (10 + 9 + 0.5 * 8) / 2.5)
        brute(L, alpha) = minimum(z + sum(max.(L .- z, 0)) / ((1 - alpha) * length(L)) for z in L)
        R = scenarios(300, 3)
        Lw = -(R * fill(1 / 3, 3))
        @test isapprox(var_cvar(Lw, 0.95)[2], brute(Lw, 0.95))
        @test_throws ArgumentError var_cvar(L, 1.0)
    end

    @testset "Euler contributions sum to the total" begin
        R = scenarios(500, 4)
        w = [0.1, 0.2, 0.3, 0.4]
        @test isapprox(sum(cvar_contributions(R, w, 0.95)), var_cvar(R, w, 0.95)[2], atol = 1e-12)
        @test isapprox(sum(vol_contributions(Sigma3, [0.2, 0.3, 0.5])), sqrt(dot([0.2, 0.3, 0.5], Sigma3, [0.2, 0.3, 0.5])))
    end

    @testset "Ledoit-Wolf" begin
        R = scenarios(250, 5)
        Sigma, delta = ledoit_wolf(R)
        @test 0 <= delta <= 1
        @test issymmetric(Sigma) && isposdef(Sigma)
        rng = MersenneTwister(7)
        big = randn(rng, 20_000, 3) .* sqrt.([1.0 2.0 3.0]) .* 0.01
        @test ledoit_wolf(big)[2] < 0.05             # lots of data: little shrinkage
    end

    @testset "min variance matches the closed form" begin
        wcf = Sigma3 \ ones(3)
        wcf ./= sum(wcf)
        @test all(>(0), wcf)                         # long-only constraint inactive
        @test isapprox(min_variance(Sigma3), wcf, atol = 1e-4)
        w = min_variance(Sigma3, Constraints(cap = 0.4))
        @test maximum(w) <= 0.4 + 1e-6 && isapprox(sum(w), 1)
    end

    @testset "tangency matches the closed form" begin
        mu, rf = [0.05, 0.09, 0.07], 0.01
        wcf = Sigma3 \ (mu .- rf)
        wcf ./= sum(wcf)
        @test all(>(0), wcf)
        @test isapprox(tangency(mu, Sigma3, rf), wcf, atol = 1e-4)
        @test_throws ArgumentError tangency(mu, Sigma3, 0.2)
        sharpe(w) = (dot(mu, w) - rf) / sqrt(dot(w, Sigma3, w))
        wc = tangency(mu, Sigma3, rf, Constraints(cap = 0.4))
        @test maximum(wc) <= 0.4 + 1e-6
        @test sharpe(wc) <= sharpe(wcf) + 1e-9
        @test sharpe(wc) >= sharpe(fill(1 / 3, 3)) - 1e-9
    end

    @testset "min CVaR linear program" begin
        R = scenarios(500, 4)
        w, obj = min_cvar(R, 0.95)
        @test isapprox(obj, var_cvar(R, w, 0.95)[2], atol = 1e-7)
        @test obj <= var_cvar(R, fill(0.25, 4), 0.95)[2] + 1e-10
        @test obj <= var_cvar(R, min_variance(ledoit_wolf(R)[1]), 0.95)[2] + 1e-10
        wc, objc = min_cvar(R, 0.95, Constraints(cap = 0.3))
        @test maximum(wc) <= 0.3 + 1e-6 && objc >= obj - 1e-10
        @test_throws ArgumentError min_cvar(R, 0.95, Constraints(cap = 0.2))
    end

    @testset "frontiers are monotone and meet their targets" begin
        R = scenarios(400, 4)
        mu = vec(mean(R, dims = 1))
        Sigma = ledoit_wolf(R)[1]
        t, W = mv_frontier(mu, Sigma, Constraints(); n = 8)
        v = [dot(w, Sigma, w) for w in eachcol(W)]
        @test all(diff(v) .>= -1e-10)
        @test all(W' * mu .>= t .- 1e-7)
        t, W = cvar_frontier(R, 0.95; n = 8)
        cv = [var_cvar(R, w, 0.95)[2] for w in eachcol(W)]
        @test all(diff(cv) .>= -1e-8)
        @test all(W' * mu .>= t .- 1e-7)
    end

    @testset "analyze" begin
        R = scenarios(600, 4)
        w = [0.4, 0.3, 0.2, 0.1]
        res = analyze(R, w; rf = 0.02, cap = 0.5)
        names = [p.name for p in res.portfolios]
        @test names[1] == "current" && names[end] == "min CVaR, mean >= current"
        cur, best = res.portfolios[1], res.portfolios[end]
        @test best.stats.cvar <= cur.stats.cvar + 1e-9
        @test dot(res.mu, best.w) >= dot(res.mu, w) - 1e-7
        @test isapprox(sum(res.contributions.cvar), cur.stats.cvar, atol = 1e-12)
        @test occursin("CVaR95", sprint(print_report, ["A", "B", "C", "D"], res))
    end

    @testset "files" begin
        dir = mktempdir()
        prices = joinpath(dir, "prices.csv")
        write(prices, "date,AAA,BBB\n2024-01-02,10,50\n2024-01-03,10.5,49\n2024-01-04,10.2,51\n")
        dates, tickers, P = read_prices(prices)
        @test tickers == ["AAA", "BBB"] && length(dates) == 3
        @test isapprox(simple_returns(P), [0.05 -0.02; 10.2 / 10.5 - 1 51 / 49 - 1])
        weights = joinpath(dir, "weights.txt")
        write(weights, "# money held\nAAA 600\nBBB 400  # euro\n")
        @test isapprox(read_weights(weights, tickers), [0.6, 0.4])
        write(weights, "CCC 1\n")
        @test_throws ArgumentError read_weights(weights, tickers)
        R = scenarios(200, 2)
        out = write_frontier(joinpath(dir, "f.csv"), ["AAA", "BBB"], R; n = 5)
        lines = readlines(out)
        @test lines[1] == "frontier,mean,vol,cvar,AAA,BBB" && length(lines) == 11
    end
end
