struct SolverResults
    weights::Vector{Float64}
    converged::Bool
    msg::String
end

function callsolver(cov, b, max_iter, tol, bounds, solver::Symbol = :ccd)
    if solver == :newton
        return newton(cov, b, max_iter, tol, bounds)
    elseif solver == :ccd
        return ccd(cov, b, max_iter, tol, bounds)
    elseif solver == :fastccd
        return fastccd(cov, b, max_iter, tol, bounds)
    else
        return fastnewton(cov, b, tol, bounds)
    end
end

"""
    minimumvariance(cov, w, [max_iter], [tol], [bounds], [solver])

```julia
cov::AbstractMatrix      covariance matrix
w::AbstractVector{Float} portfolio weights
max_iter::Int=10000      number of iterations 
tol::Float=10(-4)        the minimum tolerance of the result
bounds::Bool=true        whether to run bounds checks or not
solver::Symbol=:newton   which solver to use
```

The risk budgets are set equal to the initial weights

```math
    b_i = w_i, i \\in [1, \\ldots, N].
```

External links
* Bruder, Benjamin and Roncalli, Thierry,
  Managing Risk Exposures Using the Risk Budgeting Approach (January 20, 2012),
  doi: [10.2139/ssrn.2009778 ](http://dx.doi.org/10.2139/ssrn.2009778)

* Thomas Raffinot,
  The Journal of Portfolio Management Multi-Asset Special Issue 2018, 44 (2) 89-99,
  doi: [10.3905/jpm.2018.44.2.089 ](https://doi.org/10.3905/jpm.2018.44.2.089)
"""
function minimumvariance(cov::AbstractMatrix, w::AbstractVector; max_iter::Int64 = 10000,
    tol::Float64 = 10^(-4), bounds::Bool = true, solver::Symbol=:newton)
    # Add check that length of w is num of rows of cov
    @assert all(w.>0) == true


    b = w
    return callsolver(cov, b, max_iter, tol, bounds, solver)
end

"""
    mostdiversified(cov, w, [max_iter], [tol], [bounds], [solver])

```julia
cov::AbstractMatrix      covariance matrix
w::AbstractVector{Float} portfolio weights
max_iter::Int=10000      number of iterations 
tol::Float=10(-4)        the minimum tolerance of the result
bounds::Bool=true        whether to run bounds checks or not
solver::Symbol=:newton   which solver to use
```
Given the weights set by the user, the most diversified portfolio risk budgets are proportional to the weights and the volatility as following

```math
    b_i = \\frac{w_i \\sigma_i}{\\sum_{i=1}^N w_i \\sigma_i} , \\ i \\in [1, \\ldots, N].
```

External links
* Bruder, Benjamin and Roncalli, Thierry,
  Managing Risk Exposures Using the Risk Budgeting Approach (January 20, 2012),
  doi: [10.2139/ssrn.2009778 ](http://dx.doi.org/10.2139/ssrn.2009778)

* Thomas Raffinot,
  The Journal of Portfolio Management Multi-Asset Special Issue 2018, 44 (2) 89-99,
  doi: [10.3905/jpm.2018.44.2.089 ](https://doi.org/10.3905/jpm.2018.44.2.089)
"""
function mostdiversified(cov::AbstractMatrix, w::AbstractVector; max_iter::Int64 = 10000,
    tol::Float64 = 10^(-4), bounds::Bool = true, solver::Symbol=:newton)
    # The weight vector must be positive
    @assert all(w.>0) == true
    # The covariance matrix must be NxN
    @assert (size(cov)[1] == size(cov)[2]) == true
    # Add check that length of w is num of rows of cov

    N = size(cov)[1]
    σ = sqrt.(diag(cov))
    b = (σ .* w) / (σ' * w)



    return callsolver(cov, b, max_iter, tol, bounds, solver)
end

"""
    equalriskcontribution(cov, [max_iter], [tol], [bounds], [solver])

```julia
cov::AbstractMatrix      covariance matrix
w::AbstractVector{Float} portfolio weights
max_iter::Int=10000      number of iterations 
tol::Float=10(-4)        the minimum tolerance of the result
bounds::Bool=true        whether to run bounds checks or not
solver::Symbol=:newton   which solver to use
```
The equal risk contribution portfolio sets the risk budget for each asset in the
portfolio to the inverse of the total number of assets

```math
    b_i = \\frac{1}{N}, i \\in [1, \\ldots, N].
```

External links
* Bruder, Benjamin and Roncalli, Thierry,
  Managing Risk Exposures Using the Risk Budgeting Approach (January 20, 2012),
  doi: [10.2139/ssrn.2009778 ](http://dx.doi.org/10.2139/ssrn.2009778)

* Thomas Raffinot,
  The Journal of Portfolio Management Multi-Asset Special Issue 2018, 44 (2) 89-99,
  doi: [10.3905/jpm.2018.44.2.089 ](https://doi.org/10.3905/jpm.2018.44.2.089)
"""
function equalriskcontribution(cov::AbstractMatrix; max_iter::Int64 = 10000,
    tol::Float64 = 10^(-4), bounds::Bool = true, solver::Symbol=:newton)
    # The covariance matrix must be NxN
    @assert (size(cov)[1] == size(cov)[2]) == true

    N = size(cov)[1]
    b = [1/N for i=1:N]

    return callsolver(cov, b, max_iter, tol, bounds, solver)
end

"""
    inversevariance(cov, [max_iter], [tol], [bounds], [solver])

```julia
cov::AbstractMatrix      covariance matrix
max_iter::Int=10000      number of iterations 
tol::Float=10(-4)        the minimum tolerance of the result
bounds::Bool=true        whether to run bounds checks or not
solver::Symbol=:newton   which solver to use
```
The inverse variance portfolio is defined using the inverse variance of the financial time series
```math
    b_i = \\frac{\\sigma_i^{-2}}{\\sum_{i=1}^N \\sigma_i^{-2}}, i \\in [1, \\ldots, N].
```

External links
* Bruder, Benjamin and Roncalli, Thierry,
  Managing Risk Exposures Using the Risk Budgeting Approach (January 20, 2012),
  doi: [10.2139/ssrn.2009778 ](http://dx.doi.org/10.2139/ssrn.2009778)

* Thomas Raffinot,
  The Journal of Portfolio Management Multi-Asset Special Issue 2018, 44 (2) 89-99,
  doi: [10.3905/jpm.2018.44.2.089 ](https://doi.org/10.3905/jpm.2018.44.2.089)
"""
function inversevariance(cov::AbstractMatrix; max_iter::Int64 = 10000,
    tol::Float64 = 10^(-4), bounds::Bool = true, solver::Symbol=:newton)
    # The covariance matrix must be NxN
    @assert (size(cov)[1] == size(cov)[2]) == true

    σ = sqrt.(diag(cov))
    b = σ.^(-2) ./ sum(σ.^(-2))

    return callsolver(cov, b, max_iter, tol, bounds, solver)
end
