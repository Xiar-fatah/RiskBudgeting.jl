function _covtocorr(cov::AbstractMatrix)
    σ = diagm(sqrt.(diag(cov)))
    cov_inv = inv(σ)
    return cov_inv * (cov * cov_inv)
end