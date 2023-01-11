using ReTest

module tester
    using ReTest
    include("./theta_from_cylindrical.jl")

    function manual_theta_from_cylindrical(r, amplitudes; tol = 1e-5)
        order = length(amplitudes);
        ζ(θ::Float64) = @. $sum(amplitudes * (collectPl(cos(θ), lmax = order).parent)[2:end]);
        r_from_θ(θ)   = @. sin(θ) * (1 + ζ(θ)) - r;
        for theta = pi:(-tol):pi/2
            if abs(r_from_θ(theta)) < tol || r_from_θ(theta + tol) < 0 < r_from_θ(theta - tol) 
                return theta
            end
        end
        throw("Couldn't converge!: $r, $amplitudes")
    end
    function test_amplitudes(range, amplitudes)
            for ii in range
                @test theta_from_cylindrical(ii, amplitudes) ≈ manual_theta_from_cylindrical(ii, amplitudes) atol=1e-3
            end
    end
    
    @testset "No amplitudes" begin
        @test theta_from_cylindrical(0, [0]) ≈ pi rtol=1e-4
        test_amplitudes(range(0.0, 0.99, length=15), [0, 0]);
    end

    @testset "Amplitudes 0, 0.1" begin
        amp = [0.0, 0.1];
        test_amplitudes(range(0.0, 0.7, length = 5), amp);
    end

    @testset "Amplitudes 0, 0, 0.1" begin
        amp = [0.0, 0.0, 0.1];
        test_amplitudes(range(0.0, 0.95, length = 5), amp);
    end

    @testset "Other amplitudes" begin
        amp = [ 0.0, 0.1, 0.2, 0.3];
        test_amplitudes(range(0.0, 0.95, length = 5), amp);
        amp[2] = 0.5;
        test_amplitudes(range(0.0, 0.87, length = 5), amp);
    end
end

retest();