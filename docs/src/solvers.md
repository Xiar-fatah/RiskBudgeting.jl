# [Solvers](@id solvers)
The system of equations presented in the portfolio section can be solved with either cyclical coordinate descent as suggested by Griveau-Billion etal. or by Newton's method proposed by Florin Spinu. An comparision between the methods are made in the R language by Farah Bouzida which concludes that Newton's method is faster and more robust. Furthermore, Choi, J., & Chen applies modifications to both the cyclical coordinate descent and Newton's method with the aim of making them faster. 

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
