@testset "equalriskcontribution tests" begin
    covtest = [0.3 0.1
                0.1 0.8]
    wtest = [1/2 for i = 1:size(covtest, 1)]
    @test equalriskcontribution(covtest; solver = :fastccd).weights ≈ fastccd(covtest, wtest).weights
    @test equalriskcontribution(covtest; solver = :ccd).weights ≈ ccd(covtest, wtest).weights
    @test equalriskcontribution(covtest; solver = :newton).weights ≈ newton(covtest, wtest).weights

end

@testset "minimumvariance tests" begin
    covtest = [0.3 0.1
                0.1 0.8]

    wminus = [-1/2 for i = 1:size(covtest, 1)]
    @test_throws AssertionError minimumvariance(covtest, wminus)


end

@testset "mostdiversified tests" begin
    covtest = [0.3 0.1
                0.1 0.8]

    wminus = [-1/2 for i = 1:size(covtest, 1)]
    @test_throws AssertionError mostdiversified(covtest, wminus)


end
