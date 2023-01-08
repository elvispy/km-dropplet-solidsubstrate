using ReTest

module PSF_tester
    using ReTest, ClassicalOrthogonalPolynomials, LinearAlgebra, Polynomials
    include("./problemStructsAndFunctions.jl")
    using .problemStructsAndFunctions: LP, sinPoly;
    
    s = -1:0.02:1; NMAX = 150;
    X = LP(NMAX);
    @testset "Test Legendre Polynomial Generation" begin
        for d = 1:15:(NMAX-1) 
            @test norm(X[d].(s) - legendrep.(d, s)) ≈ 0 atol=1e-7 # 
        end
    end

    @testset "Testing sinPoly" begin
        @test norm(sin.(s) - sinPoly.(s)) ≈ 0 atol=1e-8
    end

    @testset "Testing Antiderivatives" begin
        # https://www.wolframalpha.com/input?i=int_%7Bpi%2F3%7D%5E%7Bpi%7D+sin%28x%29+*+x+*+x+dx
        Y = integrate(sinPoly * X[1]^2);
        @test (Y(pi) - Y(pi/3)) ≈ (-3 -pi/sqrt(3) + 19 * pi^2 / 18) rtol=1e-4

        # https://www.wolframalpha.com/input?i=int_%7B-1%7D%5E%7Bcos%28pi%2F3%29%7D+x+*+1%2F2+*+%283x%5E2-1%29+dx
        Y2 = integrate(X[1] * X[2]); 
        @test (Y2(cos(pi/3)) - Y2(-1)) ≈ -21/128 rtol=1e-6

        # https://www.wolframalpha.com/input?i=int_-1%5E%7Bcos%280.5%29%7D+1%2F2+*+%283x%5E2+-+1%29+*+1%2F128+%28315+x+-+4620+x%5E3+%2B+18018+x%5E5+-+25740+x%5E7+%2B+12155+x%5E9%29+dx
        Y3 = integrate(X[2] * X[9]);
        @test (Y3(cos(0.5)) - Y3(-1)) ≈ 0.00855822 atol =1e-5
        
        # Legendre Polynomials are Orthogonal  
        Y4 = integrate(X[49] * X[119]);
        @test Y4(-1) ≈ Y4(1) rtol=1e-5

    end
end

retest()