# [Solvers](@id solvers)
The system of equations presented in the portfolio section can be solved with either cyclical coordinate descent as suggested by Griveau-Billion etal. or by Newton's method proposed by [Florin Spinu](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=2297383). An comparision between the methods are made in the R language by Farah Bouzida which concludes that Newton's method is faster and more robust. Furthermore, Choi, J., & Chen applies modifications to both the cyclical coordinate descent and Newton's method with the aim of making them faster. In this package you can call the solver indepedently as shown in this example

```jldoctest basics
julia> using RiskBudgeting

julia> cov = [0.1 0.3
        0.3 0.8];

julia> b = [1/2 for i= 1:2];

julia> retval = ccd(cov, b);

julia> print(round.(retval; digits=4))
[0.7388, 0.2612]





```

## [Cyclical coordinate descent](@id solver_ccd)


```@docs
ccd
```

```@docs
fastccd
```

## [Newton's method](@id solver_newton)

```@docs
newton
```

!!! note
    Note that the modified Spinu's algorithm (fastnewton) proposed by Choi, J., & Chen, R. can obtain negative weights, which according to the constraints of risk budgeting should not occur. In addition, the identical
    solver the authors use, the hybrj routine, does not have a stable version in Julia at this moment.

```@docs
fastnewton
```
