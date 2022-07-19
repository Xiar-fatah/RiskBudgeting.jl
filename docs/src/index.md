# RiskBudgeting
*RiskBudgeting.jl* is a Julia package for calculating the weights of the risk budgeting portfolio. It consists of two steps, the risk partition of the portfolio assets and the solver for the optimization of the risk partition.


## Getting started
To begin using RiskBudgeting.jl, install it by typing
```julia
using Pkg; Pkg.add("https://github.com/Xiar-fatah/RiskBudgeting.jl")
```
Riskbudgeting.jl consists of two partitions:
- [Portfolios](@ref portfolios_intro), ready to use risk budgeting portfolios from Thomas Raffinot's Hierarchical Clustering-Based Asset Allocation. This includes equal risk contribution, inverse variance and most diversified portfolios. You then have the opportunity to call the desired solver for that portfolio or use the default setting.

- [CCD](@ref solver_ccd) and [Newton's Method](@id solver_newton) contains as of now, four implementations variations of solvers in reference to Farah Bouzida's Robust Risk Budgeting Algorithms in R and Jaehyuk Choi and Rong Chen Improved iterative methods for solving risk parity portfolio.

## Citing
Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.
