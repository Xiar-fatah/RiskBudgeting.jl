@testset "ccd tests" begin
    b_test = [0.25 for i=1:3]
    cov_test = [1 0.5 0.5
        0.5 1 0.5
        0.5 0.5 1]

    @test ccd(cov_test, b_test) ≈ [1/3 for i=1:3]

    b_test = [1/3 for i=1:3]
    @test ccd(cov_test, b_test) ≈ [1/3 for i=1:3]

    b_minus = [-1/3 for i = 1:3]
    @test_throws AssertionError ccd(cov_test, b_minus)

    cov_non_square = [1 0.5
                0.5 1
                0.5 0.5]
    @test_throws AssertionError ccd(cov_non_square, b_test)




end

@testset "fastccd tests" begin
    b_test = [0.25 for i=1:3]
    cov_test = [1 0.5 0.5
        0.5 1 0.5
        0.5 0.5 1]

    @test isapprox(fastccd(cov_test, b_test), [1/3 for i=1:3]; atol = 0.1)

    b_test = [1/3 for i=1:3]
    #@test fastccd(cov_test, b_test) ≈ [1/3 for i=1:3]

    b_minus = [-1/3 for i = 1:3]
    @test_throws AssertionError fastccd(cov_test, b_minus)

    cov_non_square = [1 0.5
                0.5 1
                0.5 0.5]
    @test_throws AssertionError fastccd(cov_non_square, b_test)

end
@testset "newton tests" begin
    b_test = [0.25 for i=1:3]
    cov_test = [1 0.5 0.5
        0.5 1 0.5
        0.5 0.5 1]

    #@test fastccd(cov_test, b_test) ≈ [1/3 for i=1:3]

    b_test = [1/3 for i=1:3]
    #@test fastccd(cov_test, b_test) ≈ [1/3 for i=1:3]

    b_minus = [-1/3 for i = 1:3]
    @test_throws AssertionError newton(cov_test, b_minus)

    cov_non_square = [1 0.5
                0.5 1
                0.5 0.5]
    @test_throws AssertionError newton(cov_non_square, b_test)




end

@testset "fastnewton tests" begin
    b_test = [0.25 for i=1:3]
    cov_test = [1 0.5 0.5
        0.5 1 0.5
        0.5 0.5 1]

    #@test fastccd(cov_test, b_test) ≈ [1/3 for i=1:3]

    b_test = [1/3 for i=1:3]
    #@test fastccd(cov_test, b_test) ≈ [1/3 for i=1:3]

    b_minus = [-1/3 for i = 1:3]
    @test_throws AssertionError newton(cov_test, b_minus)

    cov_non_square = [1 0.5
                0.5 1
                0.5 0.5]
    @test_throws AssertionError newton(cov_non_square, b_test)


end