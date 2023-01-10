using LegendrePolynomials

include("./theta_from_cylindrical.jl");
include("./problemConditionStruct.jl");
# TODO: Properly document this function
"""
    getNextStep
Tries to advance one step (Δt) in the solution of the dropplet / solid substrate interaction
"""
function getNextStep(previous_conditions::Union{ProblemConditions, Vector{ProblemConditions}}, 
    new_number_contact_points::Int64, Δt::Float64, Δr::Float64, spatial_tol::Float64, PROBLEM_CONSTANTS
    )::Tuple{ProblemConditions, Float64}

    if previous_conditions <: ProblemConditions; previous_conditions = [previous_conditions]; end

    # TODO: calculate maximum number of contact points
    if new_number_contact_points < 0 || new_number_contact_points > 100 
        errortan = Inf;
        probable_next_conditions = ProblemConditions(NaN, 
            [NaN], [NaN], [NaN], NaN, NaN, NaN, NaN, NaN
        )
    else
        # Try with same pressure distribution
        probable_next_conditions = ProblemConditions(NaN, [NaN], [NaN], 
            previous_conditions[end].pressure_amplitudes, NaN, NaN, NaN, NaN, NaN);
        #pressure_amplitudes_tentative = current_conditions.pressure_amplitudes;
        iteration = 0;
        previous_tentatives = [];
        while iteration < 100
            iteration = iteration + 1;

            # Try to advance time
            probable_next_conditions = advance_conditions(probable_next_conditions, 
                previous_conditions, new_number_contact_points, Δt, PROBLEM_CONSTANTS);
            
            probable_next_conditions, is_it_acceptable, previous_tentatives = 
                update_tentative_heuristic(probable_next_conditions, previous_conditions, 
                    Δr, Δt, spatial_tol, PROBLEM_CONSTANTS; previous_tentatives = previous_tentatives);
            
            if is_it_acceptable == true
                break;
            end
            if iteration == 100; println("Hey!") end
        end

        errortan = NaN; # TODO: Update this

    end # End main if

    return probable_next_conditions, errortan

end # end main function definition


function advance_conditions(probable_next_conditions::ProblemConditions, previous_conditions::Vector{ProblemConditions},  
    new_nb_contact_points::Int64, Δt::Float64, PROBLEM_CONSTANTS
    )::ProblemConditions

    # if previous_conditions <: ProblemConditions; previous_conditions = [previous_conditions]; end
    n = length(previous_conditions); # Determines the order of the method
    1 ≤ n ≤ 2 ? nothing : throw("Only implicit euler and BDF-2 were implemented");

    nb_harmonics = current_conditions.nb_harmonics
    pressure_amplitudes_tentative = probable_next_conditions.pressure_amplitudes;

    # Deformation amplitudes in the new coordinate system at all necessary times
    Y  = zeros(Float64, n+1, nb_harmonics, 2); 
    # Independent term (pressure term) in new coordinates and at all necessary times. (nb of times used, nb of harmonics, Spatial dimension)
    PB = zeros(Float64, 1, nb_harmonics, 2); 

    amplitudes_tent = zeros(Float64, (nb_harmonics, ));
    amplitudes_velocities_tent = zeros(Float64, (nb_harmonics, ));

    # Array of coefs sutch that coeffs .* (exp(tD) * Y) ≈ Δt × PB^{k+1}
    if n == 1
        coefs = [-1.0, 1.0];
    elseif n == 2
        rk = Δt/previous_conditions[end].dt;
        ak = (1+2rk)/(1+rk);
        bk = -(1+rk);
        ck = rk^2/(1+rk);
        coefs = [ck, bk, ak]; 
    end

    for ii = 2:nb_harmonics # Do I really need to translate into new coordinates?
        #Translating the problem into new coordinates
        PB[1, ii, :] = PROBLEM_CONSTANTS["ODE_inverse_matrices"][:, :, ii] * 
            [0; -ii * pressure_amplitudes_tentative[ii]];
        
        for jj = 1:n
            Y[jj, ii, :] = PROBLEM_CONSTANTS["ODE_inverse_matrices"][:, :, ii] * 
                [previous_conditions[jj].deformation_amplitudes[ii]; previous_conditions[jj].velocities_amplitudes[ii]];
        end

        # Extracting tentative solution
        ω_i = PROBLEM_CONSTANTS["omegas_frequencies"][ii];
        diagexp(dt::Float64) = exp(diagm(-[-1.0im * dt * ω_i, 1.0im * dt * ω_i]));
        dt(jj::Int) = Δt + previous_conditions[end].current_time - previous_conditions[jj].current_time;

        rhs_vector = zeros(Float64, (2, ));
        for jj = 1:n
            rhs_vector = rhs_vector + coefs[jj] * diagexp(dt(jj)) * Y[jj, ii, :]
        end
        Y[end, ii, :] = (Δt * PB[1, ii, :]  - rhs_vector)/ coefs[end] 
         
        Y[end, ii, :] = PROBLEM_CONSTANTS["ODE_matrices"][:, :, ii] * Y[end, ii, :];        
    end

    amplitudes_tent = Y[end, :, 1];
    amplitudes_velocities_tent = Y[end, :, 2];

    exct(jj::Int, field::String) = previous_conditions[jj].Symbol(field);
    new_CM_velocity_times_dt = - sum(coefs[1:n] .* exct.(1:n, "center_of_mass_velocity")) * Δt;

    Cl(l::Integer) = nb_harmonics >= l >= 1 ? l*(l-1)/(2*l-1) : 0;
    Dl(l::Integer) = nb_harmonics >= l >= 1 ? (l+2) * (l+1) / (2*l + 3) : 0;

    new_CM_velocity_times_dt +=  - Δt^2 * PROBLEM_CONSTANTS["froude_nb"]; 
    #Special case: First harmonics
    new_CM_velocity_times_dt +=  - Δt^2 * pressure_amplitudes_tentative[1] 
    # General case
    for ll = 2:(nb_harmonics-1)
        new_CM_velocity_times_dt += 3 * Δt^2 *
            amplitudes_tent[ll] / (2*ll+1) * (
                Cl(ll) * pressure_amplitudes_tentative[ll-1] - 
                Dl(ll) * pressure_amplitudes_tentative[ll+1]
            ) 
    end

    # Special case: Last harmonic
    new_CM_velocity_times_dt += 3 * Δt^2 *  
        amplitudes_tent[nb_harmonics] / (2 * nb_harmonics + 1) * (
            Cl(nb_harmonics) * pressure_amplitudes_tentative[nb_harmonics - 1]
        ) 

    new_CM_velocity_times_dt /= coefs[end];
    
    new_center_of_mass = (new_CM_velocity_times_dt - sum(coefs[1:n] .* exct.(1:n, "center_of_mass")));
    new_center_of_mass /= coefs[end];

    return ProblemConditions(
        nb_harmonics,
        amplitudes_tent,
        amplitudes_velocities_tent,
        pressure_amplitudes_tentative,
        current_conditions.current_time + Δt,
        Δt,
        new_center_of_mass,
        new_CM_velocity_times_dt/Δt, # Divide by Δt!
        new_nb_contact_points
    )
end


"""
    update_tentative(::ProblemConditions, ::Float64)::Tuple{ProblemConditions, Bool}

Returns a probable next solution and a flag boolean to decide whether the new solution is acceptable or not
"""
function update_tentative(probable_next_conditions::ProblemConditions, previous_conditions::Vector{ProblemConditions}, 
    Δr::Float64, Δt::Float64, spatial_tol::Float64, PROBLEM_CONSTANTS::Dict
    )::Tuple{ProblemConditions, Bool}

    harmonics_qtt = probable_next_conditions.nb_harmonics;
    is_it_acceptable = true
    # First, lets check if the given probable next condition is acceptable or not.
    heights = fill(probable_next_conditions.center_of_mass, (probable_next_conditions.number_contact_points, ));
   
    for ii = 1:probable_next_conditions.number_contact_points
        # Angle of last contact point
        θ = theta_from_cylindrical(Δr*(ii-1), probable_next_conditions.deformation_amplitudes)
        heights[ii] += cos(θ) * (1 + sum(amplitudes .* (collectPl(cos(θ), lmax = order).parent))); # We need .parent as collectPl returns an offset zero-indexed offsetarray!
        if abs(heights[ii]) > spatial_tol
            is_it_acceptable = false;
            break;
        end
    end
    if is_it_acceptable == false
        # Tentative one: Modify the discrete amplitudes, find out the new amplitudes and the consequent pressure distribution.
        θ = theta_from_cylindrical(Δr * (probable_next_conditions.number_contact_points - 1), 
            probable_next_conditions.deformation_amplitudes);

        int_endpoint = cos(θ);
        
        #d = 1 + sum([probable_next_conditions.deformation_amplitudes[ii] * (mod(ii, 2) == 0 ? 1 : -1)
             #   for ii = 1:nb_harmonics]);
        intermediate_amplitudes = zeros(Float64, (harmonics_qtt, ));
        for ii = 2:harmonics_qtt
            # Integrating from cos(θ) to 1 (known amplitude)
            intermediate_amplitudes[ii] = PROBLEM_CONSTANTS["LEGENDRE_POLYNOMIALS"][ii](1) - PROBLEM_CONSTANTS["LEGENDRE_POLYNOMIALS"][ii](int_endpoint)
            
            new_tentative_amplitude_1 = probable_next_conditions.deformation_amplitudes .* # TODO CHANGE THIS 
                 [(PROBLEM_CONSTANTS["poly_antiderivatives"][ii, jj](1) 
                    - PROBLEM_CONSTANTS["poly_antiderivatives"][ii, jj](int_endpoint)) # ∫P_ii × P_jj dx  
                        for jj = 2:harmonics_qtt] 
            intermediate_amplitudes[ii] += sum(new_tentative_amplitude_1) 

            # Integrating from -1 to cos(θ) (flat surface)
            g = PROBLEM_CONSTANTS["LPdX"][ii] # ∫ P_ii/x dx
            intermediate_amplitudes[ii] += -probable_next_conditions.center_of_mass * (g(int_endpoint) - g(-1)); 
        end

        # Now lets extract the pressure distribution that would give these intermediate amplitudes
        # Array of coefs sutch that coeffs .* (exp(tD) * Y) ≈ Δt × PB^{k+1}
        if n == 1
            coefs = [-1.0, 1.0];
        elseif n == 2
            rk = Δt/previous_conditions[end].dt;
            ak = (1+2rk)/(1+rk);
            bk = -(1+rk);
            ck = rk^2/(1+rk);
            coefs = [ck, bk, ak]; 
        end

    end

    return probable_next_conditions, is_it_acceptable
end


"""
    update_tentative_heuristic(::ProblemConditions,  ::Vector{ProblemConditions}, 
    ::Float64, ::Float64, ::Float64, ::Dict)::Tuple{ProblemConditions, Bool}

Returns a probable next solution and a flag boolean to decide whether the new solution is acceptable or not
Using an heuristic
"""
function update_tentative_heuristic(probable_next_conditions::ProblemConditions, previous_conditions::Vector{ProblemConditions}, 
    Δr::Float64, Δt::Float64, spatial_tol::Float64, PROBLEM_CONSTANTS::Dict; previous_tentatives::Vector{Vector{Float64}}
    )::Tuple{ProblemConditions, Bool}

    harmonics_qtt = probable_next_conditions.nb_harmonics;
    is_it_acceptable = true
    # First, lets check if the given probable next condition is acceptable or not.
    heights = fill(probable_next_conditions.center_of_mass, (probable_next_conditions.number_contact_points, ));
   
    for ii = 1:probable_next_conditions.number_contact_points
        # Angle of last contact point
        θ = theta_from_cylindrical(Δr*(ii-1), probable_next_conditions.deformation_amplitudes)
        heights[ii] += cos(θ) * (1 + sum(amplitudes .* (collectPl(cos(θ), lmax = order).parent))); # We need .parent as collectPl returns an offset zero-indexed offsetarray!
        if abs(heights[ii]) > spatial_tol
            is_it_acceptable = false;
            break;
        end
    end
    if is_it_acceptable == false
        # Heuristic tentative: Increase of reduce the pressure at given points to fit flat area.
        θ = theta_from_cylindrical(Δr * (probable_next_conditions.number_contact_points - 1), 
            probable_next_conditions.deformation_amplitudes);
            
        y_velocity(θ::Float64)   = sin(θ) * cos(θ) * sum(probable_next_conditions.velocities_amplitudes .* 
                (collectdnPl(cos(θ); lmax = order, n = 1).parent));
        #amp(θ::Float64) = 1 + ζ(θ);
        #r_positions = Δr * (0:(probable_next_conditions.new_number_contact_points-1));
        #position_samples = [amp(theta_from_cylindrical(r)) for r in r_positions];
        if length(previous_tentatives) == 0
            # Modify pressure amplitudes so as to try to flatten the surface

        elseif length(previous_tentatives) == 1

        else

        end

    end

    return probable_next_conditions, is_it_acceptable, previous_tentatives
end


function advance_conditions_dep(probable_next_conditions::ProblemConditions, current_conditions::ProblemConditions,  
    new_nb_contact_points::Int64, Δt::Float64, PROBLEM_CONSTANTS
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
        amplitudes_velocities_tent[ii] = Y_tent[2, ii];
    end

    # TODO: Make this integration exact assuming Al and Bl linear.
    Δt = current_conditions.dt;

    new_CM_velocity = current_conditions.center_of_mass_velocity - 
        Δt / PROBLEM_CONSTANTS["froude_nb"] - 
        Δt * (pressure_amplitudes_tentative[1] + current_conditions.pressure_amplitudes[1]) / 2;

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
