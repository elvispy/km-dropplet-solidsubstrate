using ReTest

module tester
    using ReTest
    include("./maximum_contact_radius.jl")

    function manual_max_radius(amplitudes; tol = 1e-5)
        # using bissection method
        order = length(amplitudes);
        ζ(θ::Float64)::Float64 = sum(amplitudes .* (collectPl(cos(θ), lmax = order))[2:end])
        drdθ(θ::Float64)::Float64 = cos(θ) * (1 + ζ(θ)) - 
            sin(θ)^2 * sum(amplitudes .* (collectdnPl(cos(θ), lmax = order, n = 1))[2:end]);
        
        n = 1;
        θ_min = pi/2-pi/8;
        θ_max = pi;
        @assert drdθ(θ_min) < 0 "Bad guess on θ_min!"
        @assert drdθ(θ_max) ≥ 0 "Bad guess on θ_max!"
        while b-a > tol && n < 200
            θ_mid = (θ_min+θ_max)/2;
            if abs(drdθ(θ_mid)) < tol
                return r_from_spherical(θ_mid, amplitudes);
            elseif drdθ(θ_mid) > 0
                θ_max = θ_mid;
            elseif drdθ(θ_mid) < 0
                θ_min = θ_mid;
            end
             n += 1;
        end
        throw("Couldn't converge!: $r, $amplitudes")
    end
    # println(manual_max_radius)
    function test_amplitudes(range, amplitudes)
            for ii in range
                @test theta_from_cylindrical(ii, amplitudes) ≈ manual_theta_from_cylindrical(ii, amplitudes) atol=1e-3
            end
    end
    
    @testset "No amplitudes" begin
        @test maximum_contact_radius([0.0, 0.0]) ≈ 1 rtol=1e-4
        # test_amplitudes(range(0.0, 0.99, length=15), [0, 0]);
    end
    """
    @testset "Amplitudes 0, 0.1" begin
        amp = [0.0, 0.1];
        test_amplitudes(range(0.0, 0.7, length = 5), amp);
    end
    """

end

retest();