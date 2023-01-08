using LegendrePolynomials

include("./theta_from_cylindrical.jl");
# TODO: Properly document this function
"""
    getNextStep
Tries to advance one step (Δt) in the solution of the dropplet / solid substrate interaction
"""
function getNextStep(current_conditions::ProblemConditions, new_number_contact_points::Int64, Δt::Float64, 
    Δr::Float64, spatial_tol::Float64, PROBLEM_CONSTANTS
    )::Tuple{ProblemConditions, Float64}

    if new_number_contact_points < 0 || new_number_contact_points > 100 # TODO: 
        errortan = Inf;
        probable_next_conditions = ProblemConditions(NaN, 
            [NaN], [NaN], [NaN], NaN, NaN, NaN, NaN, NaN
        )
    else
        # Try with same pressure distribution
        probable_next_conditions = ProblemConditions(NaN, [NaN], [NaN], 
            current_conditions.pressure_amplitudes, NaN, NaN, NaN, NaN, NaN);
        #pressure_amplitudes_tentative = current_conditions.pressure_amplitudes;
        iteration = 0;

        while iteration < 100
            iteration = iteration + 1;

            # Try to advance time
            probable_next_conditions = advance_conditions(current_conditions, 
                                                        PROBLEM_CONSTANTS, probable_next_conditions, new_number_contact_points, Δt);
            
            probable_next_conditions, is_it_acceptable = update_tentative(probable_next_conditions, Δr, Δt, spatial_tol, PROBLEM_CONSTANTS);
            
            if is_it_acceptable == true
                break;
            end
            if iteration == 100; println("Hey!") end
        end

        errortan = NaN; # TODO: Update this

    end # End main if

    return probable_next_conditions, errortan

end # end main function definition


function advance_conditions(current_conditions::ProblemConditions, PROBLEM_CONSTANTS, 
    probable_next_conditions::ProblemConditions, new_nb_contact_points::Int64, Δt::Float64
    )::ProblemConditions


    pressure_amplitudes_tentative = probable_next_conditions.pressure_amplitudes;

    nb_harmonics = current_conditions.nb_harmonics
    Y_old = zeros(2, nb_harmonics); # Deformation amplitudes in the new coordinate system but at time t.
    C_tent = zeros(2, nb_harmonics); # Independent term (pressure term) in new coordinates and at time t + Δt.
    C_old = zeros(2, nb_harmonics);  # Independent term in new coordinates and at time t
    Y_tent = zeros(2, nb_harmonics); # Deformation amplitudes in the new coordinate system but at time t + Δt.
    amplitudes_tent = zeros(1, nb_harmonics);
    amplitudes_velocities_tent = zeros(1, nb_harmonics);

    for ii = 2:nb_harmonics

        # Translating problem variables to new coordinates

        Y_old[:, ii] = PROBLEM_CONSTANTS["ODE_inverse_matrices"][:, :, ii] * 
            [current_conditions.deformation_amplitudes[ii]; current_conditions.velocities_amplitudes[ii]];
        C_tent[:, ii] = PROBLEM_CONSTANTS["ODE_inverse_matrices"][:, :, ii] * 
            [0; -ii * pressure_amplitudes_tentative[ii]];
        C_old[:, ii]  = PROBLEM_CONSTANTS["ODE_inverse_matrices"][:, :, ii] * 
            [0; -ii * current_conditions.pressure_amplitudes[ii]];

        # Extracting tentative solution
        ω_i = PROBLEM_CONSTANTS["omegas_frequencies"][ii];
        Y_tent[:, ii] = exp(diagm(-[-1.0im * Δt * ω_i, 1.0im * Δt * ω_i])) * 
                (Y_old + C_old * Δt / 2) + C_tent * Δt / 2;

        Y_tent[:, ii] = PROBLEM_CONSTANTS["ODE_matrices"][:, :, ii] * Y_tent[:, ii];
        
        amplitudes_tent[ii] = Y_tent[1, ii];
        amplitudes_velocities_tent = Y_tent[2, ii];
    end

    # TODO: Make this integration exact assuming Al and Bl linear.
    Δt = current_conditions.dt;

    new_CM_velocity = current_conditions.center_of_mass_velocity - 
        Δt / PROBLEM_CONSTANTS["froude_nb"] - 
        Δt * (pressure_amplitudes_tentative[1] + current_conditions.pressure_amplitudes[1]) / 2;
;
    #new_center_of_mass = current_conditions.center_of_mass + 
    #    Δt * current_conditions.center_of_mass_velocity - 
    #    Δt^2 / PROBLEM_CONSTANTS["froude_nb"] - 
    #    Δt^2 * (pressure_amplitudes_tentative[1] + current_conditions.pressure_amplitudes[1]) / 2;
    
    Cl(l::Integer) = nb_harmonics >= l >= 1 ? l*(l-1)/(2*l-1) : 0;
    Dl(l::Integer) = nb_harmonics >= l >= 1 ? (l+2) * (l+1) / (2*l + 3) : 0;

    for ii = 2:(nb_harmonics-1)
        new_CM_velocity =  new_CM_velocity  + 4 * pi * Δt *
            ( 
                amplitudes_tent[ii] / (2*ii+1) * (
                    Cl(ii) * pressure_amplitudes_tentative[ii-1] - 
                    Dl(ii) * pressure_amplitudes_tentative[ii+1]
                    ) + 
                current_conditions.deformation_amplitudes[ii] / (2*ii+1) * (
                    Cl(ii) * current_conditions.pressure_amplitudes[ii-1] - 
                    Dl(ii) * current_conditions.pressure_amplitudes[ii+1]
                )
            )/2
    end

    # Special case: Last harmonic
    new_CM_velocity = new_CM_velocity + Δt * 4 * pi * (
        amplitudes_tent[nb_harmonics] / (2 * nb_harmonics + 1) * (
            Cl(nb_harmonics) * pressure_amplitudes_tentative[nb_harmonics - 1]
        ) + 
        current_conditions.deformation_amplitudes[nb_harmonics] / (2 * nb_harmonics + 1) * (
            Cl(nb_harmonics) * current_conditions.pressure_amplitudes[nb_harmonics - 1]
        )
    )

    new_center_of_mass = current_conditions.center_of_mass + Δt * new_CM_velocity;

    return ProblemConditions(
        nb_harmonics,
        amplitudes_tent,
        amplitudes_velocities_tent,
        pressure_amplitudes_tentative,
        current_conditions.current_time + Δt,
        Δt,
        new_center_of_mass,
        new_CM_velocity,
        new_nb_contact_points
    )
end


"""
    update_tentative(::ProblemConditions, ::Float64)::Tuple{ProblemConditions, Bool}

Returns a probable next solution and a flag to decide wether the new solution is acceptable or not
"""
function update_tentative(probable_next_conditions::ProblemConditions, Δr::Float64, 
    Δt::Float64, spatial_tol::Float64, PROBLEM_CONSTANTS::Dict
    )::Tuple{ProblemConditions, Bool}

    harmonics_qtt = probable_next_conditions.nb_harmonics;
    is_it_acceptable = true
    heights = fill(probable_next_conditions.center_of_mass, (probable_next_conditions.number_contact_points, ));

    for ii = 1:probable_next_conditions.number_contact_points
        θ = theta_from_cylindrical(Δr*(ii-1), probable_next_conditions.deformation_amplitudes)
        heights[ii] = heights[ii] +  cos(θ) * (sum(amplitudes .* (collectPl(cos(θ), lmax = order).parent)) + 1);
        if abs(heights[ii]) > spatial_tol
            is_it_acceptable = false;
            break;
        end
    end
    if is_it_acceptable == false
        # Tentative one: Modify the discrete amplitudes, find out the new amplitudes and the consequent pressure distribution.
        θ = theta_from_cylindrical(Δr * probable_next_conditions.number_contact_points, 
                                        probable_next_conditions.deformation_amplitudes);
        
        int_endpoint = cos(θ);

        intermediate_amplitudes = zeros(Float64, (harmonics_qtt, ));
        for ii = 2:harmonics_qtt

            # Integrating from -1 (pi) to cos(θ) (θ)
            intermediate_amplitudes[ii] = PROBLEM_CONSTANTS["LEGENDRE_POLYNOMIALS"][ii+1](int_endpoint) - PROBLEM_CONSTANTS["LEGENDRE_POLYNOMIALS"][ii+1](-1)
            intermediate_amplitudes[ii] = sum([PROBLEM_CONSTANTS["poly_antiderivatives"][ii, jj] for jj = 4:10]) # TODO FIX THIS

        end
    end

    return probable_next_conditions, is_it_acceptable




end