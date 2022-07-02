"""
    ccd(cov, b, [max_iter], [tol])

    ```julia
    cov::AbstractMatrix Covariance matrix
    b::AbstractVector{Float} Risk budgeting vector
    max_iter::Int Number of iterations for cyclical coordinate descent
    tol::Float The minimum tolerance of the result
    ```
    Calculates the solution of the risk budgeting portfolio with Cyclical Coordinate Descent given the covariance matrix
    and the risk partitions between the assets.

    External links
    * Griveau-Billion, Théophile and Richard, Jean-Charles and Roncalli, Thierry,
    A Fast Algorithm for Computing High-Dimensional Risk Parity Portfolios (September 1, 2013). 
    doi: http://dx.doi.org/10.2139/ssrn.2325255 
"""
function ccd(cov::AbstractMatrix, b::AbstractVector{Float64}, 
    max_iter::Int64 = 10000, tol::Float64 = 10^(-4))::AbstractVector
    
    # The risk budgeting vector must be positive
    @assert all(b.>0) == true
    # The covariance matrix must be NxN
    @assert (size(cov)[1] == size(cov)[2]) == true

    x = 1 ./ diag(cov)
    Σx = cov * x
    xΣx = x' * Σx
    for iter = 1:max_iter
        i = (iter) % size(cov)[1] + 1
        nonsqrt = -Σx[i] + x[i]*cov[i,i]
        x̃ = (nonsqrt + sqrt(nonsqrt^2 + 4*cov[i,i] * b[i] * sqrt(xΣx))) / (2*cov[i, i])
        
        Σx = Σx + cov[:,i]*(x̃ - x[i]) # Σx̃
        xΣx = xΣx + cov[i,i] * (x[i]^2 - x̃^2) - 2*x[i]*(cov[i, :]' * x)  # σ(x̃)
        x[i] = x̃

        xΣx = xΣx + 2*x[i]* (cov[i,:]' * x)
        if maximum(abs.(x .* (cov * x) .- b)) < tol
            return x
        end
    end

    println("Cyclical Coordinate Descent has failed to converge!")
    return x  
    
end

