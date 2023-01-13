#using ReTest

# module tester_gnt
    using ReTest
    include("./get_next_step.jl")
    include("./problemStructsAndFunctions.jl");
    using .problemStructsAndFunctions: LP
    include("./problemConditionStruct.jl");
    using Plots

    # First, testing calculate_exit_angle
    function manual_exit_angle(amplitudes, θ; ϵ = 1e-9)
        order = length(amplitudes);
        ζ = zeta(amplitudes; order = order);
        # ζ(θ::Float64) = @. $sum(amplitudes * (collectPl(cos(θ), lmax = order).parent)[2:end]);
        r_from_θ(θ) = sin(θ) * (1 + ζ(θ));
        z_from_θ(θ) = cos(θ) * (1 + ζ(θ));
        return (z_from_θ(θ - ϵ) - z_from_θ(θ + ϵ)) / (r_from_θ(θ - ϵ) - r_from_θ(θ + ϵ))
    end
 
    @testset "Testing exit angle function: No amplitudes" begin
        @test calculate_exit_angle([0.0], round(pi; digits = 20)) ≈ 0 atol=1e-6
        for angle in LinRange(pi/1.9, pi/1.1, 14)
            @test calculate_exit_angle([0.0], angle) ≈ manual_exit_angle([0.0], angle) rtol = 1e-5
            @test testing_perpendicular([0.0], angle) ≈ 0 atol = 1e-8;     
        end
    end

    @testset "Testing exit angle: Amplitudes  [0, 0.2]" begin
        amp = [0.0, 0.2]

        for angle in LinRange(pi/1.9, pi, 15)
            @test calculate_exit_angle(amp, angle) ≈ manual_exit_angle(amp, angle) atol=1e-3
        end
    end

    @testset "Testing exit angle: Amplitudes  [0.2, 0.1]" begin
        amp = [0.2, 0.1]

        for angle in LinRange(pi/1.9, pi/1.1, 15)
            @test calculate_exit_angle(amp, angle) ≈ manual_exit_angle(amp, angle) rtol=1e-5
        end
    end

    @testset "Testing exit angle: Amplitudes  0.1.*[.1, .1, .1, .1, .1, .1]" begin
        amp = 0.1.*[.1, .1, .1, .1, .1, .1]

        for angle in LinRange(pi/1.9, pi/1.1, 15)
            @test calculate_exit_angle(amp, angle) ≈ manual_exit_angle(amp, angle) rtol=1e-5
        end
    end

    @testset "Testing exit angle: Amplitudes 0.15 .* [.2, -.1, .4, -.3]" begin
        amp = 0.15 .* [.2, -.1, .4, -.3];
        for angle in LinRange(pi/1.9, pi/1.1, 15)
            @test calculate_exit_angle(amp, angle) ≈ manual_exit_angle(amp, angle) rtol=1e-3
        end
    end 

    ## Now, testing advance_conditions
    # Trying with Bl = 0, for all l and times.
    harmonics_qtt = 5;
    f(n::Integer) = sqrt(n .* (n+2) .* (n-1));
    omegas_frequencies = f.(1:harmonics_qtt)';
    LEGENDRE_POLYNOMIALS = LP(harmonics_qtt);
    Δt = 0.01;

    ODE_matrices = ones(ComplexF64, 2, 2, harmonics_qtt); # Y' = -PDP^-1 Y + B ==> (exp(tD)*Y)' = e^(tD) P^-1 B;
    # Where D = diag([-im*ω_i, im*ω_i])
    ODE_matrices[2, 1, :] =  1.0im * omegas_frequencies;
    ODE_matrices[2, 2, :] = -1.0im * omegas_frequencies;

    ODE_inverse_matrices = 1/2 * ones(ComplexF64, 2, 2, harmonics_qtt);
    ODE_inverse_matrices[1, 2, :] = -0.5im ./ omegas_frequencies;
    ODE_inverse_matrices[2, 2, :] =  0.5im ./ omegas_frequencies;

    # println(ODE_matrices[:, :, 3] * [[-omegas_frequencies[3] *1.0im, 0.0];; [0.0, 1.0im * omegas_frequencies[3]]] * ODE_inverse_matrices[:, :, 3])

    PROBLEM_CONSTANTS = Dict(
        "froude_nb" => 1,#"weber_nb"  => weber_nb, #"ball_mass" => mS/mass_unit,
        "omegas_frequencies" => (f.(1:harmonics_qtt)).^2,
        "ODE_matrices" => ODE_matrices,
        "ODE_inverse_matrices" => ODE_inverse_matrices, #"poly_antiderivatives" => polynomials_antiderivatives,#"LEGENDRE_POLYNOMIALS" => LEGENDRE_POLYNOMIALS,#"LPdX" => LPdX
    )

    previous_conditions = [
        ProblemConditions(
            harmonics_qtt,
            zeros(Float64, (harmonics_qtt, )),
            0.1 *ones(Float64, (harmonics_qtt, )),
            zeros(Float64, (harmonics_qtt, )),
            0.0,
            Δt,
            0.0,
            0.0,
            0
        )
    ]

    probable_next_conditions = deepcopy(previous_conditions[1]);
    N = 100;
    CM = zeros(Float64, (N, ));
    amps = zeros(Float64, (N, harmonics_qtt));
    for ii = 1:N
        global probable_next_conditions
        probable_next_conditions = advance_conditions(probable_next_conditions, previous_conditions,
        0, Δt, PROBLEM_CONSTANTS);
        CM[ii] = probable_next_conditions.center_of_mass;
        amps[ii, :] = probable_next_conditions.deformation_amplitudes;
        
        @assert norm(probable_next_conditions.pressure_amplitudes) ≈ 0 "Pressures are not zero!";
    end

    x = Δt .* (0:(N-1));

    CM_exact = map((t) -> -t^2, x);
    amps_exact = map((t) -> 
        0.1 * sin(PROBLEM_CONSTANTS["omegas_frequencies"][3] * t), x);
    plt = plot(x, [CM, CM_exact], label = ["Numeric" "Exact"]);

    display(plt)

    #readline()
    



# end

# retest();