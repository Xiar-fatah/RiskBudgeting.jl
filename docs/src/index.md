# RiskBudgeting
*RiskBudgeting.jl* is a Julia package for calculating the weights of the risk budgeting portfolios, with more emphasis on the solvers. Riskbudgeting.jl consists of two partitions:

- [Portfolios](@ref portfolios_intro), ready to use risk budgeting portfolios from Thomas Raffinot's Hierarchical Clustering-Based Asset Allocation. This includes equal risk contribution, inverse variance and most diversified portfolios. You then have the opportunity to call the desired solver for that portfolio or use the default setting.

- [CCD](@ref solver_ccd) and [Newton's Method](@ref solver_newton) contains as of now, four implementations variations of solvers in reference to Farah Bouzida's Robust Risk Budgeting Algorithms in R and Jaehyuk Choi and Rong Chen Improved iterative methods for solving risk parity portfolio.


## Getting started
To begin using RiskBudgeting.jl, install it by typing
```julia
using Pkg; Pkg.add("RiskBudgeting")
```
or from Julia's package manager by typing `] add RiskBudgeting`.



## Citing
If you are using this package for research and wish cite this package you can do this by following this [link](https://zenodo.org/badge/latestdoi/504891771).
