# [Portfolios](@id portfolios_intro)
In comparison to the previous methods, risk budgeting portfolios considers the amount of risk contribution each asset in the portfolio can hold.  Let $w_i$ and ``\\mathcal{R}(w_1, \\hdots, w_N)`` denote the weights of the asset ``i`` and the risk measure of the portfolio, for example volatility, respectively. Given that the risk measure is coherent and convex it will satisfy the condition of Euler decomposition

```math
\\mathcal{R}(w_1, \\hdots, w_N) = \\sum_{i=1}^n w_i \\frac{\\partial \\mathcal{R}(w_1 , \\ldots, w_n)}{\\partial w_i}.
```

As a result, the risk contribution, interpreted as the amount of risk exposure an asset in a portfolio is allowed to collect, of asset $i$ is expressed as
```math
\textrm{RC}_i (w_1, \hdots, w_N) = w_i \frac{\partial \mathcal{R}(w_1 , \ldots, w_N)}{\partial w_i}.
```
This implies that if given the risk budgets $\{b_1, \ldots, b_N\}$, the risk budgeting portfolio is solved by finding the solutions to the non-linear system of equations ``b_i = \\textrm{RC}_i``. By generalizing the system of equations with the additional constraints of being long only and complete capital allocation, the risk budgeting portfolio develop into
```math
\begin{align}\label{risk_budgeting_portfolio}
    \begin{cases}
    w_i \\frac{\\partial \\mathcal{R}(w_1 , \\ldots, w_n)}{\\partial w_i} = b_i \\mathcal{R}, \\
    b_i \\geq 0, \\
    w_i \\geq 0, \\
    \\sum_{i=1}^N b_i = 1, \\
    \\sum_{i=1}^N w_i = 1.
    \\end{cases}
 \end{align}
```


## Equal risk contribution
```@docs
equalriskcontribution
```

## Inverse variance
```@docs
inversevariance
```


## Most diversified

```@docs
mostdiversified
```
