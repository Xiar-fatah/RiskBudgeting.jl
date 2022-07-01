

function ccd(cov::AbstractMatrix, b::AbstractVector, max_iter::Int64, tol::Float64)::AbstractVector
    x = 1 ./ diag(cov)
    Σx = cov * x
    xΣx = x' * Σx
    for iter = 1:max_iter
        i = iter % cov.shape[1]
        nonsqrt = -Σx[i] + x[i]*cov[i,i]
        x̃ = (nonsqrt + sqrt(nonsqrt^2 + 4*cov[i,i] * b[i] * sqrt(xΣx))) / (2*cov[i, i])
        Σx = Σx + cov[:,i]*(x_new - x[i]) # Σx̃
        xΣx = xΣx + cov[i,i] * (x[i]^2 - x̃^2) - 2*x[i]*(cov[i, :].* x)  # σ(x̃)
        x[i] = x̃
        xΣx = xΣx + 2*x[i]* (cov[i,:] .* x)
        if max(abs(x * (cov .* x) - b)) < tol
            return x
        end
    end
    println("Cyclical Coordinate Descent has failed to converge!")
    return x    
end
