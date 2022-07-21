import RiskBudgeting: ccd, fastccd, newton

function helper(cov, b, max_iter, tol, bounds, solver::Symbol)::AbstractVector
    if solver == newton
        weights = newton(cov, b, max_iter, tol, bounds)
    elseif solver == ccd
        weights = fastccd(cov, b, max_iter, tol, bounds)
    elseif solver == fastccd
        weights = ccd(cov, b, max_iter, tol, bounds)
    else
        weights = fastnewton(cov, b, max_iter, tol, bounds)
    end
    return weights
end

"""
    minimumvariance(cov, b, [max_iter], [tol], [bounds], [solver])
TBA
"""
function minimumvariance(cov::AbstractMatrix, max_iter::Int64 = 10000,
    tol::Float64 = 10^(-4), bounds::Bool = true, solver::Symbol=:newton)::AbstractVector

    return helper(cov, b, max_iter, tol, bounds, solver)
end
"""
    mostdiversified(cov, b, [max_iter], [tol], [bounds], [solver])
Given the weights set by the user, the most diversified portfolio risk budgets are proportional to the weights and the volatility as following
```math
    b_i = \\frac{w_i \\sigma_i}{\\sum_{i=1}^N w_i \\sigma_i} , \\ i \\in [1, \\ldots, N].
```

"""
# double check this, it is suppose to be set by the user.'

function mostdiversified(cov::AbstractMatrix, max_iter::Int64 = 10000,
    tol::Float64 = 10^(-4), bounds::Bool = true, solver::Symbol=:newton)::AbstractVector
    # The covariance matrix must be NxN
    @assert (size(cov)[1] == size(cov)[2]) == true

    N = size(cov)[1]
    σ = sqrt.(diag(cov))
    b = [1/N for i=1:N] 
    b = (σ .* b) / (σ' * b)

    # The risk budgeting vector must be positive
    @assert all(b.>0) == true

    return helper(cov, b, max_iter, tol, bounds, solver)
end

"""
    equalriskcontribution(cov, b, [max_iter], [tol], [bounds], [solver])

The equal risk contribution portfolio sets the risk budget for each asset in the
portfolio to the inverse of the total number of assets

```math
    b_i = \\frac{1}{N}, i \\in [1, \\ldots, N].
```
"""
function equalriskcontribution(cov::AbstractMatrix, max_iter::Int64 = 10000,
    tol::Float64 = 10^(-4), bounds::Bool = true, solver::Symbol=:newton)::AbstractVector
    # The covariance matrix must be NxN
    @assert (size(cov)[1] == size(cov)[2]) == true

    N = size(cov)[1]
    b = [1/N for i=1:N]

    # The risk budgeting vector must be positive
    @assert all(b.>0) == true
    return helper(cov, b, max_iter, tol, bounds, solver)
end

"""
    inversevariance(cov, b, [max_iter], [tol], [bounds], [solver])

The inverse variance portfolio is defined using the inverse variance of the financial time series
```math
    b_i = \\frac{\\sigma_i^{-2}}{\\sum_{i=1}^N \\sigma_i^{-2}}.
```
"""
function inversevariance(cov::AbstractMatrix, max_iter::Int64 = 10000,
    tol::Float64 = 10^(-4), bounds::Bool = true, solver::Symbol=:newton)::AbstractVector
    # The covariance matrix must be NxN
    @assert (size(cov)[1] == size(cov)[2]) == true

    σ = sqrt.(diag(cov))
    b = σ.^(-2) ./ sum(σ.^(-2))

    # The risk budgeting vector must be positive
    @assert all(b.>0) == true

    return helper(cov, b, max_iter, tol, bounds, solver)
end

