using ReTest

module tester_psf
    using ReTest, ClassicalOrthogonalPolynomials, LinearAlgebra, Polynomials
    include("./problemStructsAndFunctions.jl")
    using .problemStructsAndFunctions: LP, sinPoly, integrate_poly;
    
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
    end # End testset
    a = 1;
    @testset "Testing ∫P(x)/x^m dx" begin

        # https://www.wolframalpha.com/input?i=int_1%5E2+1%2F8+%2815+x+-+70+x%5E3+%2B+63+x%5E5%29%2Fx+dx
        @test integrate_poly(X[5], 1)(2) - integrate_poly(X[5], 1)(1) ≈ 1817/60 rtol=1e-5 #"∫ P_5(x)/x dx not properly integrated!"
        
        # https://www.wolframalpha.com/input?i=int_%7B0.5%7D%5E%7B1.2%7D+1%2F8+%283+-+30+x%5E2+%2B+35+x%5E4%29%2Fx+dx
        @test integrate_poly(X[4], 1)(1.2) - integrate_poly(X[4], 1)(0.5) ≈ 0.296691 atol=1e-5 #"∫ P_4(x)/x dx not properly integrated!"
        
        # https://www.wolframalpha.com/input?i=int_%7B0.2%7D%5E%7B0.9%7D+1%2F2+%28-3+x+%2B+5+x%5E3%29%2Fx%5E2+dx
        @test integrate_poly(X[3], 2)(0.9) - integrate_poly(X[3], 2)(0.2) ≈ -1.29362 atol=1e-5 #"∫ P_4(x)/x dx not properly integrated!"
    end
end

retest()