using ReTest

module tester
    using ReTest
    include("./maximum_contact_radius.jl")
    include("./zeta.jl");
    include("./r_from_spherical.jl");

    function manual_max_radius(amplitudes; tol = 1e-5)
        # using bissection method
        order = length(amplitudes);
        # ζ(θ::Float64)::Float64 = sum(amplitudes .* (collectPl(cos(θ), lmax = order))[2:end])
        ζ = zeta(amplitudes; order = order);
        drdθ(θ::Float64)::Float64 = cos(θ) * (1 + ζ(θ)) - 
            sin(θ)^2 * sum(amplitudes .* (collectdnPl(cos(θ), lmax = order, n = 1))[1:end]);
        
        n = 1;
        θ_min = pi/2-pi/8;
        θ_max = 3.14159;
        while drdθ(θ_min) < 0 
            θ_min = rand() * pi/2;
        end
        @assert drdθ(θ_max) ≤ 0 "Bad guess on θ_max!"
        while θ_max-θ_min > tol && n < 200
            θ_mid = (θ_min+θ_max)/2;
            if abs(drdθ(θ_mid)) < tol
                return r_from_spherical(θ_mid, amplitudes);
            elseif drdθ(θ_mid) > 0
                θ_min = θ_mid;
            elseif drdθ(θ_mid) < 0
                θ_max = θ_mid;
            end
             n += 1;
        end
        throw("Couldn't converge!: $amplitudes, $θ_max")
    end
    # println(manual_max_radius)
    function test_amplitudes(range, amplitudes)
            for ii in range
                @test theta_from_cylindrical(ii, amplitudes) ≈ manual_theta_from_cylindrical(ii, amplitudes) atol=1e-3
            end
    end
    
    @testset "No amplitudes" begin
        @test maximum_contact_radius([0.0, 0.0]) ≈ 1 rtol=1e-5
        # test_amplitudes(range(0.0, 0.99, length=15), [0, 0]);
    end
    
    @testset "Amplitudes .1, -.1, .4, .3" begin
        amp = [.1, -.1, .4, .3];
        @test maximum_contact_radius(amp) ≈ manual_max_radius(amp) atol=1e-4
    end

    @testset "Amplitudes  0.15 * (.2, -.1, .4, -.3)" begin
        amp = 0.15 .* [.2, -.1, .4, -.3];
        @test maximum_contact_radius(amp) ≈ manual_max_radius(amp) atol=1e-4
    end
    

end

retest();