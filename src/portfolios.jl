import RiskBudgeting: ccd



"""
    mostdiversified(cov, b, [max_iter], [tol])
"""
function mostdiversified(cov::AbstractMatrix, 
    max_iter::Int64 = 10000, tol::Float64 = 10^(-4))::AbstractVector
    # The covariance matrix must be NxN
    @assert (size(cov)[1] == size(cov)[2]) == true

    N = size(cov)[1]
    σ = sqrt.(diag(cov))
    b = [1/N for i=1:N]
    b = (σ .* b) / (σ' * b)

    # The risk budgeting vector must be positive
    @assert all(b.>0) == true

    return ccd(cov, b, max_iter, tol, predefined = 1)
end

"""
    equalriskcontribution(cov, b, [max_iter], [tol])
"""
function equalriskcontribution(cov::AbstractMatrix, 
    max_iter::Int64 = 10000, tol::Float64 = 10^(-4))::AbstractVector
    # The covariance matrix must be NxN
    @assert (size(cov)[1] == size(cov)[2]) == true

    N = size(cov)[1]
    b = [1/N for i=1:N]

    # The risk budgeting vector must be positive
    @assert all(b.>0) == true
    return ccd(cov, b, max_iter, tol)
end

"""
    inversevariance(cov, b, [max_iter], [tol])
"""
function inversevariance(cov::AbstractMatrix, 
    max_iter::Int64 = 10000, tol::Float64 = 10^(-4))::AbstractVector
    # The covariance matrix must be NxN
    @assert (size(cov)[1] == size(cov)[2]) == true

    σ = sqrt.(diag(cov))
    b = σ.^(-2) ./ sum(σ.^(-2))

    # The risk budgeting vector must be positive
    @assert all(b.>0) == true

    return ccd(cov, b, max_iter, tol, predefined = 1)
end