using RiskBudgeting
using Test
using LinearAlgebra


tests = ["ccd"]
println("Runing tests:")
for t in tests
    fp = "$(t).jl"
    println("* $fp ...")
    include(fp)
end

