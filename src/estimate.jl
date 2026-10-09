"""
    ledoit_wolf(R) -> (Sigma, delta)

Ledoit-Wolf (2004) shrinkage of the sample covariance towards the
scaled identity `m I`, with `m` the average sample variance:

    Sigma = delta m I + (1 - delta) S

The intensity `delta in [0, 1]` is the estimated optimal one. The sample
covariance `S` uses the `1/T` normalization, as in the paper.
"""
function ledoit_wolf(R::AbstractMatrix{<:Real})
    T, N = size(R)
    T > 1 || throw(ArgumentError("need at least two observations"))
    X = R .- mean(R, dims = 1)
    S = (X' * X) ./ T
    m = tr(S) / N
    d2 = sum(abs2, S - m * I) / N
    d2 > 0 || return (Matrix(S), 0.0)    # S is already a multiple of I
    b2bar = sum(sum(abs2, x * x' - S) for x in eachrow(X)) / (N * T^2)
    delta = min(b2bar, d2) / d2
    Sigma = delta * m * I + (1 - delta) * S
    return Matrix(Symmetric(Sigma)), delta
end
