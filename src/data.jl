"""
    read_prices(path) -> (dates, tickers, P)

Read a CSV with header `date,T1,...,Tn` and one row per date in
ascending order. Prices must be positive and complete: the merge done
by `scripts/fetch.sh` keeps only the dates common to every ticker.
"""
function read_prices(path::AbstractString)
    data, header = readdlm(path, ',', header = true)
    size(data, 2) >= 2 || throw(ArgumentError("$path: expected date plus at least one ticker"))
    tickers = String.(strip.(string.(vec(header)[2:end])))
    dates = String.(string.(data[:, 1]))
    issorted(dates) || throw(ArgumentError("$path: dates must be ascending"))
    body = data[:, 2:end]
    all(x -> x isa Real, body) || throw(ArgumentError("$path: missing or non-numeric prices"))
    P = Float64.(body)
    all(>(0), P) || throw(ArgumentError("$path: prices must be positive"))
    return dates, tickers, P
end

"""
    read_weights(path, tickers) -> w

Read lines `TICKER amount` (`#` starts a comment). Amounts can be
weights or money held: they are normalized to sum to one. Tickers not
listed get weight zero.
"""
function read_weights(path::AbstractString, tickers::AbstractVector{<:AbstractString})
    w = zeros(length(tickers))
    for (lineno, raw) in enumerate(eachline(path))
        line = strip(first(split(raw, '#'; limit = 2)))
        isempty(line) && continue
        fields = split(line)
        length(fields) == 2 || throw(ArgumentError("$path:$lineno: expected `TICKER amount`"))
        i = findfirst(==(fields[1]), tickers)
        i === nothing && throw(ArgumentError("$path:$lineno: unknown ticker $(fields[1])"))
        w[i] = parse(Float64, fields[2])
    end
    all(>=(0), w) || throw(ArgumentError("$path: short positions are not supported"))
    s = sum(w)
    s > 0 || throw(ArgumentError("$path: amounts sum to zero"))
    return w ./ s
end

"""
    simple_returns(P) -> R

Per-period simple returns, one row per period and one column per asset.
"""
simple_returns(P::AbstractMatrix{<:Real}) = P[2:end, :] ./ P[1:end-1, :] .- 1
