# Command-line front end: julia --project=. bin/cvarport.jl PRICES WEIGHTS [options]

using CVaRPort

const USAGE = """
usage: julia --project=. bin/cvarport.jl PRICES.csv WEIGHTS.txt [options]

  --alpha A      CVaR confidence level (default 0.95)
  --rf R         annual risk-free rate, e.g. 0.025 (default 0)
  --cap C        maximum weight per asset (default 1)
  --periods P    return periods per year: 252 daily, 52 weekly, 12 monthly (default 252)
  --no-shrink    sample covariance instead of Ledoit-Wolf
  --frontier F   also write both efficient frontiers to the CSV file F
"""

function parse_args(args)
    pos = String[]
    opt = Dict{String,String}()
    flags = Set{String}()
    i = 1
    while i <= length(args)
        a = args[i]
        if a in ("-h", "--help")
            print(USAGE)
            exit(0)
        elseif a == "--no-shrink"
            push!(flags, a)
        elseif a in ("--alpha", "--rf", "--cap", "--periods", "--frontier")
            i < length(args) || error("$a needs a value")
            opt[a] = args[i += 1]
        elseif startswith(a, "-")
            error("unknown option $a")
        else
            push!(pos, a)
        end
        i += 1
    end
    length(pos) == 2 || error("expected PRICES and WEIGHTS")
    return pos, opt, flags
end

function main(args)
    pos, opt, flags = parse_args(args)
    alpha = parse(Float64, get(opt, "--alpha", "0.95"))
    rf = parse(Float64, get(opt, "--rf", "0"))
    cap = parse(Float64, get(opt, "--cap", "1"))
    periods = parse(Int, get(opt, "--periods", "252"))
    shrink = !("--no-shrink" in flags)

    _, tickers, P = read_prices(pos[1])
    w = read_weights(pos[2], tickers)
    R = simple_returns(P)
    res = analyze(R, w; alpha = alpha, rf = rf, periods = periods, cap = cap, shrink = shrink)
    print_report(tickers, res)

    if haskey(opt, "--frontier")
        write_frontier(opt["--frontier"], tickers, R;
                       alpha = alpha, periods = periods, cap = cap, shrink = shrink)
        println("\nfrontiers written to ", opt["--frontier"])
    end
end

try
    main(ARGS)
catch err
    err isa Union{ArgumentError,ErrorException,DimensionMismatch,SystemError} || rethrow()
    println(stderr, "cvarport: ", sprint(showerror, err))
    occursin("expected PRICES", sprint(showerror, err)) && print(stderr, USAGE)
    exit(1)
end
