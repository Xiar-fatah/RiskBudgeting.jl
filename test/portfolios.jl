@testset "portfolio tests" begin
    cov_test = [1 0.5 0.5
        0.5 1 0.5
        0.5 0.5 1]

    println(mostdiversified(cov_test))
    #@test_throws AssertionError mostdiversified(cov_test, b_minus)

    #@test_throws AssertionError equalriskcontribution(cov_test, b_minus)


    #ArgumentError([msg])

end