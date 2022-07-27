"""
    ccd(cov, b, [max_iter], [tol], [bounds])

```julia
cov::AbstractMatrix Covariance matrix
b::AbstractVector{Float} Risk budgeting vector
max_iter::Int Number of iterations for cyclical coordinate descent
tol::Float The minimum tolerance of the result
bounds::Bool Whether to run bounds checks or not
```
Calculates the solution of the risk budgeting portfolio with cyclical coordinate descent given the covariance matrix
and the risk partitions between the assets.

External links
* Griveau-Billion, Théophile and Richard, Jean-Charles and Roncalli, Thierry,
  A Fast Algorithm for Computing High-Dimensional Risk Parity Portfolios (September 1, 2013).
  doi: [10.2139/ssrn.2325255](http://dx.doi.org/10.2139/ssrn.2325255)
"""
function ccd(cov::AbstractMatrix, b::AbstractVector{Float64},
    max_iter::Int64 = 10000, tol::Float64 = 10^(-4), bounds::Bool = true)

    if bounds == true
        # The risk budgeting vector must be positive
        @assert all(b.>0) == true
        # The covariance matrix must be NxN
        @assert (size(cov,1) == size(cov,2)) == true
    end

    σ = sqrt.(diag(cov))
    x = (ones(size(cov,1)) ./ (σ)) / (1 / sum(σ))
    x̃ = x
    Σx = cov * x
    xΣx = x' * Σx
    for iter = 1:max_iter
        i = (iter) % size(cov,1) + 1
        bb = -Σx[i] + x[i]*cov[i,i]
        x̃[i] = (bb + sqrt(bb^2 + 4*cov[i,i] * b[i] * sqrt(xΣx))) / (2*cov[i, i])
        Σx = Σx + cov[:,i]*(x̃[i] - x[i]) # Σx̃
        xΣx = xΣx + cov[i,i] * (x[i]^2 - x̃[i]^2) - 2*x[i]*(cov[i, :]' * x) + 2*x̃[i] * (cov[i,:]' * x̃)  # σ(x̃)
        rc = (x̃ .* (cov * x̃)) ./ (sqrt(x̃' * (cov * x̃)))

        if maximum(abs.(rc / sum(rc) .- b)) < tol
            weights =  x̃ / sum(x̃)
            return SolverResults(weights, true, "Cyclical coordinate descent has succeeded to converge!")
        end
        x = x̃

    end
    weights = x̃ / sum(x̃)
    return SolverResults(weights, false, "Cyclical coordinate descent has failed to converge!")
end


"""
    fastccd(cov, b, [max_iter], [tol], [bounds])

```julia
cov::AbstractMatrix Covariance matrix
b::AbstractVector{Float} Risk budgeting vector
max_iter::Int Number of iterations for cyclical coordinate descent
tol::Float The minimum tolerance of the result
bounds::Bool Whether to run bounds checks or not

```
A faster solution finder of the risk budgeting portfolio with cyclical coordinate descent given the covariance matrix
and the risk partitions between the assets.

External links
* Choi, J., & Chen, R. (2022).
  Improved iterative methods for solving risk parity portfolio,
  Journal of Derivatives and Quantitative Studies, 30(2)
  doi: [10.48550/arXiv.2203.00148](https://doi.org/10.48550/arXiv.2203.00148)
"""
function fastccd(cov::AbstractMatrix, b::AbstractVector{Float64},
    max_iter::Int64 = 10000, tol::Float64 = 10^(-4), bounds::Bool = true)


    if bounds == true
        # The risk budgeting vector must be positive
        @assert all(b.>0) == true
        # The covariance matrix must be NxN
        @assert (size(cov,1) == size(cov,2)) == true
    end


    corr = _covtocorr(cov)
    σ = sqrt.(diag(cov))
    N = size(corr,1)
    onevec = ones(N)
    x = onevec ./sqrt(sqrt(sum(corr)))

    for iter = 1:max_iter
        i = (iter) % size(cov,1) + 1

        aᵢ = 0.5 * ((corr * x)[i] - x[i])
        x[i] = sqrt(aᵢ^2 + b[i]) - aᵢ
        x = x ./ (sqrt(x' * (corr * x)))
        if maximum(abs.(x .* (corr * x) .- b)) < tol
            weights = (x./σ) ./ sum((x./σ))
            return SolverResults(weights, true, "Fast cyclical coordinate descent has succeeded to converge!")
        end
    end
    println("Cyclical Coordinate Descent has failed to converge!")
    weights = (x./σ) ./ sum((x./σ))
    return SolverResults(weights, false, "Fast cyclical coordinate descent has failed to converge!")
end

"""
    newton(cov, b, [max_iter], [tol], [bounds])

```julia
cov::AbstractMatrix Covariance matrix
b::AbstractVector{Float} Risk budgeting vector
max_iter::Int Number of iterations for cyclical coordinate descent
tol::Float The minimum tolerance of the result
bounds::Bool Whether to run bounds checks or not

```
Solution finder of the risk budgeting portfolio with Newton's method given the covariance matrix
and the risk partitions between the assets (also known as Spinu's algorithm).

External links
* Spinu, Florin,
  An Algorithm for Computing Risk Parity Weights (July 30, 2013).
  doi: [10.2139/ssrn.2297383](http://dx.doi.org/10.2139/ssrn.2297383)

"""
function newton(cov::AbstractMatrix, b::AbstractVector{Float64},
    max_iter::Int64 = 10000, tol::Float64 = 10^(-4), bounds::Bool = true)

    if bounds == true
        # The risk budgeting vector must be positive
        @assert all(b.>0) == true
        # The covariance matrix must be NxN
        @assert (size(cov,1) == size(cov,2)) == true
    end

    # Initilization
    u = ones(size(cov,1))
    x = sqrt(sum(b) / (u' * (cov * u))) * u
    λₛₜₐᵣ = 0.95 * (3-sqrt(5)) / 2
    # Damped Phase
    for i = 1:max_iter
        λₖ, Δx, δₖ = iteration(cov, x, u, b)
        x = x - Δx ./ (1+δₖ)
        if λₖ > λₛₜₐᵣ
            break
        end
    end
    # Quadratic Phase
    for i = 1:max_iter
        _, Δx, _ = iteration(cov, x, u, b)
        x = x - Δx
        rc = (x .* (cov * x)) ./ (sqrt(x' * (cov * x)))
        if maximum(abs.(rc / sum(rc) .- b)) < tol
            weights = x ./ sum(x)
            return SolverResults(weights, true, "Newton's method has succeeded to converge!")
        end
    end
    weights = x ./ sum(x)
    return SolverResults(weights, false, "Newton's method has failed to converge!")
end

function iteration(cov, x, u, b)
    uₖ = cov * x - b ./ x
    Hₖ = cov + diagm(b ./ (x .* x))
    Δx = inv(Hₖ) * uₖ
    δₖ = norm(Δx ./ x, Inf)
    λₖ = sqrt(transpose(uₖ) * Δx)
    return λₖ, Δx, δₖ
end

"""
    fastnewton(cov, b, [tol], [bounds])

```julia
cov::AbstractMatrix Covariance matrix
b::AbstractVector{Float} Risk budgeting vector
tol::Float The minimum tolerance of the result
bounds::Bool Whether to run bounds checks or not

```
A faster solution finder of the risk budgeting portfolio based on Newton's method given the covariance matrix
and the risk partitions between the assets. The core concept is that the solution will converge faster if the initial
guess of the weights is closer to the solution in comparision to Spinu's algorithm.

External links
* Choi, J., & Chen, R. (2022).
    Improved iterative methods for solving risk parity portfolio,
    Journal of Derivatives and Quantitative Studies, 30(2)
    doi: [10.48550/arXiv.2203.00148](https://doi.org/10.48550/arXiv.2203.00148)
"""
function fastnewton(cov::AbstractMatrix, b::AbstractVector{Float64},
     tol::Float64 = 10^(-4), bounds::Bool = true)
    if bounds == true
        # The risk budgeting vector must be positive
        @assert all(b.>0) == true
        # The covariance matrix must be NxN
        @assert (size(cov,1) == size(cov,2)) == true
    end
    corr = _covtocorr(cov)
    σ = diag(cov)
    onevec = ones(size(cov,1))
    a = (corr * onevec - onevec) / (2* sqrt(onevec' * corr * onevec))
    x = sqrt.(a.^2 + b) - a

    w, converged = my_solve(f!, g!, x, cov, b, tol) # Original paper used corr to test

    if converged == false
        weights = w ./ σ / sum(w ./ σ)
        return SolverResults(weights, false, "Fast Newtons method did not converge!") 
    elseif any(0 .> w) == true
        weights = w ./ σ / sum(w ./ σ)
        return SolverResults(weights, false, "One or more weights is negative, fast Newtons method failed to converge!") 
    else
        weights = w ./ σ / sum(w ./ σ)
        return SolverResults(weights, true, "Fast Newton's method has succeeded to converge!")
    end
end

function f!(fvec, x, cov, b)
    temp = cov * x - b ./ x
    fvec[1] = temp[1]
    fvec[2] = temp[2]
    fvec
end

function g!(fjac, x, cov, b)
    temp = cov + diagm(b ./ (x .* x))
    fjac[1,1] = temp[1,1]
    fjac[1,2] = temp[1,2]
    fjac[2,1] = temp[2,1]
    fjac[2,2] = temp[2,2]
    fjac
end

function my_solve(f!, g!, x, cov, b, tol)
    n = size(x,1)
    maxfev = 100 * (n+1)
    retval = fsolve((fvec,x) -> f!(fvec, x, cov, b),
        (fjac,x) -> g!(fjac, x, cov, b),
        x,
        tol = tol,
        method = :hybr,
        factor = 100.0,
        mode = 1,
        maxfev = maxfev)
    return retval.x, retval.converged
end
