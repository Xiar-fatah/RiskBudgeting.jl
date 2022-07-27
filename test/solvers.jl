@testset "ccd tests" begin
    cov_test = [1 0.5 0.5
        0.5 1 0.5
        0.5 0.5 1]
    b_test = [1/3 for i=1:3]
    b_minus = [-1/3 for i = 1:3]
    @test_throws AssertionError ccd(cov_test, b_minus)

    cov_non_square = [1 0.5
                0.5 1
                0.5 0.5]
    @test_throws AssertionError ccd(cov_non_square, b_test)




end

@testset "fastccd tests" begin
    cov_test = [1 0.5 0.5
        0.5 1 0.5
        0.5 0.5 1]
    b_test = [1/3 for i=1:3]
    b_minus = [-1/3 for i = 1:3]

    @test_throws AssertionError fastccd(cov_test, b_minus)

    cov_non_square = [1 0.5
                0.5 1
                0.5 0.5]
    @test_throws AssertionError fastccd(cov_non_square, b_test)

end
@testset "newton tests" begin
    cov_test = [1 0.5 0.5
        0.5 1 0.5
        0.5 0.5 1]
    b_test = [1/3 for i=1:3]
    b_minus = [-1/3 for i = 1:3]
    @test_throws AssertionError newton(cov_test, b_minus)

    cov_non_square = [1 0.5
                0.5 1
                0.5 0.5]
    @test_throws AssertionError newton(cov_non_square, b_test)




end

@testset "fastnewton tests" begin
    cov_test = [1 0.5 0.5
        0.5 1 0.5
        0.5 0.5 1]
    b_test = [1/3 for i=1:3]
    b_minus = [-1/3 for i = 1:3]
    @test_throws AssertionError fastnewton(cov_test, b_minus)

    cov_non_square = [1 0.5
                0.5 1
                0.5 0.5]
    @test_throws AssertionError fastnewton(cov_non_square, b_test)


end

@testset "solver tests" begin
    cov_test = [ 0.404095   0.0379728   0.007343
     0.0379728  0.407645    0.00870648
     0.007343   0.00870648  0.00113161]
    b_test = [1/3 for i=1:3]

    store_ccd = ccd(cov_test, b_test).weights
    store_fastccd = fastccd(cov_test, b_test).weights
    store_newton = newton(cov_test, b_test).weights

    @test all.(isapprox(store_newton, store_ccd,atol = 0.01))
    @test all.(isapprox(store_newton, store_fastccd,atol = 0.01))
    @test all.(isapprox(store_ccd, store_fastccd,atol = 0.01))

end
