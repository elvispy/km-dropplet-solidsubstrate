using ReTest

module tester
    using ReTest
    include("./getNextStep.jl")

    function manual_exit_angle(amplitudes, θ; ϵ = 1e-5)
        order = length(amplitudes);
        ζ = zeta(amplitudes; order = order);
        # ζ(θ::Float64) = @. $sum(amplitudes * (collectPl(cos(θ), lmax = order).parent)[2:end]);
        r_from_θ(θ) = sin(θ) * (1 + ζ(θ));
        z_from_θ(θ) = cos(θ) * (1 + ζ(θ));
        return (z_from_θ(θ - ϵ) - z_from_θ(θ + ϵ)) / (r_from_θ(θ - ϵ) - r_from_θ(θ + ϵ))
    end

    function test_amplitudes(range, amplitudes)
            for ii in range
                @test calculate_exit_angle(ii, amplitudes) ≈ manual_theta_from_cylindrical(ii, amplitudes) atol=1e-3
            end
    end
    
    @testset "Testing exit angle function: No amplitudes" begin
        @test calculate_exit_angle([0.0], round(pi; digits = 20)) ≈ 0 atol=1e-6
        for angle in LinRange(pi/1.9, pi/1.1, 14)
            @test calculate_exit_angle([0.0], angle) ≈ manual_exit_angle([0.0], angle) rtol = 1e-5
        end
    end

    @testset "Testing exit angle: Amplitudes 0.15 .* [.2, -.1, .4, -.3]" begin
        amp = 0.15 .* [.2, -.1, .4, -.3];
        for angle in LinRange(pi/1.9, pi/1.1, 15)
            @test calculate_exit_angle(amp, angle) ≈ manual_exit_angle(amp, angle) rtol = 1e-5
        end
    end


end

retest();