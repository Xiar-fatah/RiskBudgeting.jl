var documenterSearchIndex = {"docs":
[{"location":"newton/#solver_newton","page":"Newton's Method","title":"Newton's Method","text":"","category":"section"},{"location":"newton/","page":"Newton's Method","title":"Newton's Method","text":"Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.","category":"page"},{"location":"newton/","page":"Newton's Method","title":"Newton's Method","text":"newton","category":"page"},{"location":"newton/#RiskBudgeting.newton","page":"Newton's Method","title":"RiskBudgeting.newton","text":"newton(cov, b, [max_iter], [tol], [bounds])\n\ncov::AbstractMatrix Covariance matrix\nb::AbstractVector{Float} Risk budgeting vector\nmax_iter::Int Number of iterations for cyclical coordinate descent\ntol::Float The minimum tolerance of the result\nbounds::Bool Whether to run bounds checks or not\n\n\nSolution finder of the risk budgeting portfolio with Newton's method given the covariance matrix and the risk partitions between the assets (also known as Spinu's algorithm).\n\nExternal links\n\nSpinu, Florin, An Algorithm for Computing Risk Parity Weights (July 30, 2013). doi: 10.2139/ssrn.2297383\n\n\n\n\n\n","category":"function"},{"location":"newton/","page":"Newton's Method","title":"Newton's Method","text":"note: Note\nNote that the modified Spinu's algorithm (fastnewton) proposed by Choi, J., & Chen, R. can obtain negative weights, which according to the constraints of risk budgeting should not occur. In addition, the identical solver the authors use, the hybrj routine, does not have a stable version in Julia at this moment.","category":"page"},{"location":"newton/","page":"Newton's Method","title":"Newton's Method","text":"fastnewton","category":"page"},{"location":"newton/#RiskBudgeting.fastnewton","page":"Newton's Method","title":"RiskBudgeting.fastnewton","text":"fastnewton(cov, b, [tol], [bounds])\n\ncov::AbstractMatrix Covariance matrix\nb::AbstractVector{Float} Risk budgeting vector\ntol::Float The minimum tolerance of the result\nbounds::Bool Whether to run bounds checks or not\n\n\nA faster solution finder of the risk budgeting portfolio based on Newton's method given the covariance matrix and the risk partitions between the assets. The core concept is that the solution will converge faster if the initial guess of the weights is closer to the solution in comparision to Spinu's algorithm.\n\nExternal links\n\nChoi, J., & Chen, R. (2022).   Improved iterative methods for solving risk parity portfolio,   Journal of Derivatives and Quantitative Studies, 30(2)   doi: 10.48550/arXiv.2203.00148\n\n\n\n\n\n","category":"function"},{"location":"ccd/#solver_ccd","page":"Cyclical Coordinate Descent","title":"Cyclical Coordinate Descent","text":"","category":"section"},{"location":"ccd/","page":"Cyclical Coordinate Descent","title":"Cyclical Coordinate Descent","text":"Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.","category":"page"},{"location":"ccd/","page":"Cyclical Coordinate Descent","title":"Cyclical Coordinate Descent","text":"ccd","category":"page"},{"location":"ccd/#RiskBudgeting.ccd","page":"Cyclical Coordinate Descent","title":"RiskBudgeting.ccd","text":"ccd(cov, b, [max_iter], [tol], [bounds])\n\ncov::AbstractMatrix Covariance matrix\nb::AbstractVector{Float} Risk budgeting vector\nmax_iter::Int Number of iterations for cyclical coordinate descent\ntol::Float The minimum tolerance of the result\nbounds::Bool Whether to run bounds checks or not\n\nCalculates the solution of the risk budgeting portfolio with Cyclical Coordinate Descent given the covariance matrix and the risk partitions between the assets.\n\nExternal links\n\nGriveau-Billion, Théophile and Richard, Jean-Charles and Roncalli, Thierry, A Fast Algorithm for Computing High-Dimensional Risk Parity Portfolios (September 1, 2013). doi: 10.2139/ssrn.2325255\n\n\n\n\n\n","category":"function"},{"location":"ccd/","page":"Cyclical Coordinate Descent","title":"Cyclical Coordinate Descent","text":"fastccd","category":"page"},{"location":"ccd/#RiskBudgeting.fastccd","page":"Cyclical Coordinate Descent","title":"RiskBudgeting.fastccd","text":"fastccd(cov, b, [max_iter], [tol], [bounds])\n\ncov::AbstractMatrix Covariance matrix\nb::AbstractVector{Float} Risk budgeting vector\nmax_iter::Int Number of iterations for cyclical coordinate descent\ntol::Float The minimum tolerance of the result\nbounds::Bool Whether to run bounds checks or not\n\n\nA faster solution finder of the risk budgeting portfolio with Cyclical Coordinate Descent given the covariance matrix and the risk partitions between the assets.\n\nExternal links\n\nChoi, J., & Chen, R. (2022). Improved iterative methods for solving risk parity portfolio, Journal of Derivatives and Quantitative Studies, 30(2) doi: 10.48550/arXiv.2203.00148\n\n\n\n\n\n","category":"function"},{"location":"examples/#example","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.","category":"page"},{"location":"portfolios/#portfolios_intro","page":"Portfolios","title":"Portfolios","text":"","category":"section"},{"location":"portfolios/","page":"Portfolios","title":"Portfolios","text":"In comparison to the previous methods, risk budgeting portfolios considers the amount of risk contribution each asset in the portfolio can hold.  Let w_i and mathcalR(w_1 hdots w_N) denote the weights of the asset i and the risk measure of the portfolio, for example volatility, respectively. Given that the risk measure is coherent and convex it will satisfy the condition of Euler decomposition","category":"page"},{"location":"portfolios/","page":"Portfolios","title":"Portfolios","text":"mathcalR(w_1 hdots w_N) = sum_i=1^n w_i fracpartial mathcalR(w_1  ldots w_n)partial w_i","category":"page"},{"location":"portfolios/","page":"Portfolios","title":"Portfolios","text":"As a result, the risk contribution, interpreted as the amount of risk exposure an asset in a portfolio is allowed to collect, of asset i is expressed as","category":"page"},{"location":"portfolios/","page":"Portfolios","title":"Portfolios","text":"textrmRC_i (w_1 hdots w_N) = w_i fracpartial mathcalR(w_1  ldots w_N)partial w_i","category":"page"},{"location":"portfolios/","page":"Portfolios","title":"Portfolios","text":"This implies that if given the risk budgets b_1 ldots b_N, the risk budgeting portfolio is solved by finding the solutions to the non-linear system of equations b_i = textrmRC_i. By generalizing the system of equations with the additional constraints of being long only and complete capital allocation, the risk budgeting portfolio develop into","category":"page"},{"location":"portfolios/","page":"Portfolios","title":"Portfolios","text":"beginalignedlabelrisk_budgeting_portfolio\n    w_i fracpartial mathcalR(w_1  ldots w_n)partial w_i = b_i mathcalR \n    b_i geq 0 \n    w_i geq 0 \n    sum_i=1^N b_i = 1 \n    sum_i=1^N w_i = 1\n endaligned","category":"page"},{"location":"portfolios/#Equal-risk-contribution","page":"Portfolios","title":"Equal risk contribution","text":"","category":"section"},{"location":"portfolios/","page":"Portfolios","title":"Portfolios","text":"equalriskcontribution","category":"page"},{"location":"portfolios/#RiskBudgeting.equalriskcontribution","page":"Portfolios","title":"RiskBudgeting.equalriskcontribution","text":"equalriskcontribution(cov, b, [max_iter], [tol], [bounds], [solver])\n\n\n\n\n\n","category":"function"},{"location":"portfolios/#Inverse-variance","page":"Portfolios","title":"Inverse variance","text":"","category":"section"},{"location":"portfolios/","page":"Portfolios","title":"Portfolios","text":"inversevariance","category":"page"},{"location":"portfolios/#RiskBudgeting.inversevariance","page":"Portfolios","title":"RiskBudgeting.inversevariance","text":"inversevariance(cov, b, [max_iter], [tol], [bounds], [solver])\n\n\n\n\n\n","category":"function"},{"location":"portfolios/#Most-diversified","page":"Portfolios","title":"Most diversified","text":"","category":"section"},{"location":"portfolios/","page":"Portfolios","title":"Portfolios","text":"mostdiversified","category":"page"},{"location":"portfolios/#RiskBudgeting.mostdiversified","page":"Portfolios","title":"RiskBudgeting.mostdiversified","text":"mostdiversified(cov, b, [max_iter], [tol], [bounds], [solver])\n\ndouble check this, it is suppose to be set by the user.\n\n\n\n\n\n","category":"function"},{"location":"#RiskBudgeting","page":"Introduction","title":"RiskBudgeting","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"RiskBudgeting.jl is a Julia package for calculating the weights of the risk budgeting portfolio. It consists of two steps, the risk partition of the portfolio assets and the solver for the optimization of the risk partition.","category":"page"},{"location":"#Getting-started","page":"Introduction","title":"Getting started","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"To begin using RiskBudgeting.jl, install it by typing","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"using Pkg; Pkg.add(\"https://github.com/Xiar-fatah/RiskBudgeting.jl\")","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Riskbudgeting.jl consists of two partitions:","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Portfolios, ready to use risk budgeting portfolios from Thomas Raffinot's Hierarchical Clustering-Based Asset Allocation. This includes equal risk contribution, inverse variance and most diversified portfolios. You then have the opportunity to call the desired solver for that portfolio or use the default setting.\nCCD and Newton's Method contains as of now, four implementations variations of solvers in reference to Farah Bouzida's Robust Risk Budgeting Algorithms in R and Jaehyuk Choi and Rong Chen Improved iterative methods for solving risk parity portfolio.","category":"page"},{"location":"#Citing","page":"Introduction","title":"Citing","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.","category":"page"}]
}
