# [Solvers](@id solvers)
Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.

## [Cyclical Coordinate Descent](@id solver_ccd)


```@docs
ccd
```

```@docs
fastccd
```

# [Newton's Method](@id solver_newton)

```@docs
newton
```

!!! note
    Note that the modified Spinu's algorithm (fastnewton) proposed by Choi, J., & Chen, R. can obtain negative weights, which according to the constraints of risk budgeting should not occur. In addition, the identical
    solver the authors use, the hybrj routine, does not have a stable version in Julia at this moment.

```@docs
fastnewton
```
