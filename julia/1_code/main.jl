"""
    Author: Elvis Aguero.
    Date: 11th December, 2022.
    main
"""

module km_main

# Exports
export solveMotion


#using Printf;
using Plots;
default(legend = false);
using JLD2, FileIO;
using Polynomials;
#using Dates;
#using CSV;
#using DataFrames;


#Imports
include("./problemStructsAndFunctions.jl");
using .problemStructsAndFunctions # Specialized structs and functions to be used in this file
#--- Nothing here


# Main function
"""
    solveMotion:
        Tries to solve the kinematic match problem between a deformable droplet and a 
        solid substrate under axistymmetric conditions.
"""

function solveMotion(; # <== Keyword Arguments! 
    # OBS: All units are in CGS!
    undisturbed_radius::Float64 = .238,  # Radius of the undeformed spherical sphere 
    initial_height::Float64 = Inf,    # Initial position of the sphere center of mass of the sphere (Inf = start barely touching)
    initial_velocity::Float64 = -1.0, # Initial velocity of the sphere 
    initial_amplitude::Vector{Float64} = Float64[], # Initial amplitudes of the dropplet (Default = undisturbed) OBS: First index is A_1
    rhoS::Float64 = NaN,            # Sphere's density
    sigmaS::Float64 = NaN,          # Sphere's Surface Tension
    g::Float6 = 9.8065e+2,          # Gravitational constant
    harmonics_qtt::Int64 = 20,      # Number of harmonics to be used 
    nb_pressure_samples::Int64 = NaN,      # Number of intervals in contact radius (NaN = Equal to number of harmonics)
    max_dt::Float64 = 0.01,         # maximum allowed temporal time step
    angle_tol::Float64 = 5/360 * 2 * pi, # Angle tolerance to accept a solution (in radians) 
    spatial_tol::Float64 = 1e-8,    # Tolerance to accept that dropplet touches the substrate
    simulation_time::Float64 = 10.0,# Maximum allowed time
    live_plotting::Bool = false     # Whether to plot or not the live results
    )

    if isnan(nb_pressure_samples)
        spatial_step = 1/harmonics_qtt;
    else
        spatial_step = 1/nb_pressure_samples;
    end

    # Dimensionless Units
    time_unit = undisturbed_radius/initial_velocity; # Temporal dimensionless number
    length_unit = rS;
    velocity_unit = abs(initial_velocity);
    pressure_unit = rhoS * velocity_unit^2;
    froude_nb   = initial_velocity.^2/(g*undisturbed_radius);
    weber_nb    = rhoS * undisturbed_radius * initial_velocity.^2 / sigmaS; # Weber's number of the dropplet
    reynolds_nb = undisturbed_radius * velocity_unit / nu; # Reynolds' number
    mS = 4*pi*undisturbed_radius^3 * rhoS / 3;
    mass_unit = rhoS * length_unit^3;

    ## Initial conditions
    # Set dropplet's sphere height initial conditions
    if initial_height == Inf
        initial_height = get_initial_height(initial_amplitude)
    else
        @assert initial_height >= get_initial_height(initial_amplitude)
        initial_height = initial_height/length_unit;
    end
    if length(initial_amplitude) == 1
        initial_amplitude = zeros((harmonics_qtt, ));
    else
        initial_amplitude = initial_amplitude/length_unit;
    end

    if length(amplitudes_velocities) == 1
        amplitudes_velocities = zeros((harmonics_qtt, ));
    else
        amplitudes_velocities = amplitudes_velocities/velocity_unit;
    end

    initial_pressure_coefficients = zeros((harmonics_qtt, )) / pressure_unit; # Just to emphasize the units of these coefficients.

    if max_dt == 0
        Δt = 0.01; 
    else 
        Δt = max_dt/time_unit; 
    end

    initial_time = 0;
    current_time = initial_time/time_unit;
    final_time = simul_time/time_unit;
    current_index = 2; # This integer points to the next available index in variables that are going to 
                       # export data (index 1 is for initial conditions)
    maximum_index = ceil(Int64, (final_time - initial_time)/Δt) + 4;
    number_of_extra_indexes = 0;

    contact_points = 0; # Initial number of contact points
    contact_time = 0; #TODO: Lab contact time and so on?
    grow_dt = false; # THis variable controls how fast Δt can grow
    iii = 0; jjj = 0; # Indexes to keep track how small is Δt compared to max_dt

    #= if isfile("./2_pipeline/LEGENDRE_POLYNOMIALS.jld2") == true
        LEGENDRE_POLYNOMIALS = load_object("LEGENDRE_POLYNOMIALS.jld2");
        n = length(LEGENDRE_POLYNOMIALS);
        if n > harmonics_qtt
            LEGENDRE_POLYNOMIALS = LEGENDRE_POLYNOMIALS[1:harmonics_qtt];
        else:
     =#   
    LEGENDRE_POLYNOMIALS = LP(harmonics_qtt);
    LPdX = Vector{Function}(undef, harmonics_qtt);
    polynomials_antiderivatives = Matrix{Function}(undef, harmonics_qtt, harmonics_qtt);

    # This function calculates the primitive of P/x:  (P::Polynomial) -> ∫P(x) * x^(-1) dx
    function integrate_poly(p::Polynomial{Float64, :x})::Function
        y = integrate(Polynomial([p.coeffs[ii] for ii = 2:degree(p)]));
        return (t) -> y(t) + p.coeffs[1] * log(t)
    end

    @variables x
    for ii = 1:harmonics_qtt
        for jj = 1:ii
            # This data structure at [ii, jj] will return the integral of P_{ii-1} * P_{jj-1}
            polynomials_antiderivatives[ii, jj] = integrate( LEGENDRE_POLYNOMIALS[ii] * LEGENDRE_POLYNOMIALS[jj])
        end
        # This array has the integral of P_{ii-1}/x
        LPdX[ii] = integrate_poly(LEGENDRE_POLYNOMIALS[ii]);
    end

    # https://www.wolframalpha.com/input?i=int_1%5E2+1%2F8+%2815+x+-+70+x%5E3+%2B+63+x%5E5%29%2Fx+dx
    @assert isapprox(LPdX[5](2) - LPdX[5](1), 1817/60; rtol=1e-5) "∫ P_5(x)/x dx not properly integrated!"
    
    # https://www.wolframalpha.com/input?i=int_%7B0.5%7D%5E%7B1.2%7D+1%2F8+%283+-+30+x%5E2+%2B+35+x%5E4%29%2Fx+dx
    @assert isapprox(LPdX[4](1.2) - LPdX[4](0.5), 0.296691; atol=1e-5) "∫ P_4(x)/x dx not properly integrated!"
    

    f(n::Integer) = sqrt(n .* (n+2) .* (n-1) ./ weber_nb);

    omegas_frequencies = f(1:harmonics_qtt)';

    ODE_matrices = zeros(2, 2, harmonics_qtt); # Y' = -PDP^-1 Y + B ==> (exp(tD)*Y)' = e^(tD) P^-1 B;
    ODE_matrices[1, 1, :] =  ones(1, harmonics_qtt);
    ODE_matrices[1, 2, :] =  ones(1, harmonics_qtt);
    ODE_matrices[2, 1, :] =  1.0im * omegas_frequencies;
    ODE_matrices[2, 2, :] = -1.0im * omegas_frequencies;

    ODE_inverse_matrices = 1/2 * ones(2, 2, N);
    ODE_inverse_matrices[1, 2, :] = -0.5im ./ omegas_frequencies;
    ODE_inverse_matrices[2, 2, :] =  0.5im ./ omegas_frequencies;

    PROBLEM_CONSTANTS = Dict(
        "froude_nb" => froude_nb,
        "weber_nb"  => weber_nb, 
        "ball_mass" => mS/mass_unit;
        "ODE_matrices" => ODE_matrices,
        "ODE_inverse_matrices" => ODE_inverse_matrices,
        "poly_antiderivatives" => polynomials_antiderivatives,
        "LEGENDRE_POLYNOMIALS" => LEGENDRE_POLYNOMIALS,
        "LPdX" => LPdX
    )

    probable_next_conditions = Vector{ProblemConditions}(undef, 5);
    current_conditions = ProblemConditions(harmonics_qtt, initial_amplitude, 
            amplitudes_velocities, initial_pressure_coefficients, 0.0, Δt, 
            initial_height, initial_velocity, 0);
    previous_conditions = [current_conditions, current_conditions]; # TODO: Define this array properly to implement BDF2.

    if ("julia" in readdir());  cd("julia\\");  end
    if ("1_code" in readdir()); cd("1_code\\"); end
    A = match(r"pipeline", file_name)
    if A !== nothing
        #A = A[1]; 
        SAVE_PATH = dirname(file_name);
    else
        SAVE_PATH = "../2_pipeline/default/out";
        file_name = "$SAVE_PATH/$file_name";
    end

    if isdir(SAVE_PATH) == false; mkpath(SAVE_PATH); end
    if (isfile(file_name) == false) && false # TODO: Decide what to export
        headers_data_frame = DataFrame(
            ID  = []
        )
        CSV.write(file_name, headers_data_frame)
    end

    if plotter == true
        plot_width = ceil(Int64, 3 * N);
        xplot = LinRange(0, plot_width/N, plot_width);
        ηX = [-xplot[end:-1:2]; xplot] * length_unit;
        # ηU  = zeros(1, 2 * plot_width - 1);
        #step = ceil(Int64, N/15);

        θ = LinRange(0, 2*π, 500)
        circleX = length_unit * sin.(θ);
        circleY = length_unit * cos.(θ);
        # For animation:
        # https://stackoverflow.com/questions/58103056/erasing-previous-data-plots-in-julia-plots-jl-gr-backend

    end

    ## Preparing post-processing

    # Preallocate variables that will be exported (All of them have units!)
    recorded_height = zeros((maximum_index, )); recorded_height[1] = initial_height * length_unit;
    recorded_pressure = zeros(maximum_index, nb_pressure_samples);
    recorded_contact_points = zeros(Int64, (maximum_index,)); recorded_contact_points[1] = 0;
    recorded_velocity = zeros((maximum_index, )); recorded_velocity[1] = initial_velocity * velocity_unit;
    recorded_times = zeros((maximum_index, )); recorded_times[1] = initial_time * time_unit;
    # TODO: implement the following version:
    recorded_conditions = Vector{ProblemConditions}(undef, (maximum_index, )); recorded_conditions[1] = previous_conditions[end];

    # Coefficient of restitution
    mechanical_energy_in = NaN;
    mechanical_energy_out = NaN; # TODO: Lab COef of restitution?


    while ( current_time < final_time)
        errortan = Inf * ones((5, ));
        recalculate = false;

        # First, we try to solve with the same number of contact points
        probable_next_conditions[3], errortan[3] = getNextStep(previous_conditions, contact_points, Δt, spatial_step, 
                spatial_tol, PROBLEM_CONSTANTS);

        if abs(errortan[3]) < 1e-8 # If almost no error, we accept the solution
            current_conditions = probable_next_conditions[3];
            previous_conditions = [previous_conditions[2:end]..., probable_next_conditions[3]];
        else # If there is some error, we try with diffferent contact points
            # Lets try with one more point
            probable_next_conditions[4], errortan[4] = getNextStep(previous_conditions, contact_points + 1, Δt, spatial_step, 
            spatial_tol, PROBLEM_CONSTANTS);
            # Lets try with one point less
            probable_next_conditions[2], errortan[2] = getNextStep(previous_conditions, contact_points - 1, Δt, spatial_step, 
            spatial_tol, PROBLEM_CONSTANTS);
              
            if (abs(errortan[3]) > abs(errortan[4]) ||  (abs(errortan[3]) > abs(errortan[2])))
                if abs(errortan[4]) <= abs(errortan[2])
                    # Now lets check with one more point to be sure
                    _, errortan[5] = getNextStep(previous_conditions, contact_points + 2, Δt, spatial_step, 
                    spatial_tol, PROBLEM_CONSTANTS);

                    if abs(errortan[4]) < abs(errortan[5])
                        # Accept new data 
                        previous_conditions = [previous_conditions[2:end]..., probable_next_conditions[4]];
                        current_conditions  = probable_next_conditions[4];
                        contact_points      = contact_points + 1;
                    else
                        recalculate = true;
                    end
                else
                    # now lets check if errortan is good enough with one point less
                    _, errortan[1] = getNextStep(previous_conditions, contact_points - 2, Δt, spatial_step, 
                    spatial_tol, PROBLEM_CONSTANTS);

                    if abs(errortan[2]) < abs(errortan[1])
                        # Accept new data
                        previous_conditions = [previous_conditions[2:end]..., probable_next_conditions[2]];
                        current_conditions  = probable_next_conditions[2];
                        contact_points      = contact_points - 1;
                    else
                        recalculate = true;
                    end 

                end # End of errortan[4] <= errortan[2]
            else # The same number of contact points may be the best
                if errortan[3] == Inf # ==> All errors are infinity
                    recalculate = true;
                else
                    # Accept new data
                    previous_conditions = [previous_conditions[2:end]..., probable_next_conditions[3]];
                    current_conditions  = probable_next_conditions[3];
                end
            end
        end # End outer if 

        if recalculate == true
            Δt = Δt/2;
            # Refine time step in index notation 
            iii = iii + 1; jjj = 2 * jjj;
        else
            current_time = current_time + Δt; jjj = jjj + 1;
            if mod(jjj, 2) == 0 && grow_dt == true
                jjj = div(jjj, 2); 
                iii = iii - 1;
                # Increase time step
                Δt = 2 * Δt;
                # Decrease the number of time you can make dt bigger
                grow_dt = false;
            end

            #  TODO: Update Indexes if necessary

            # TODO: # Stored data

            # If we are in a multiple of max_dt, reset indexes
            if jjj == 2^iii
                jjj = 0;
                grow_dt = true;
                #indexes_to_save[current_to_save] = current_index - 1;
                #current_to_save = current_to_save + 1;
            else
                number_of_extra_indexes = number_of_extra_indexes + 1;
            end

            if live_plotting == true
                # Do some live plotting here
            else
                # Do some real-time variable updating here
            end

        end


    end # end main while loop


end # end main function declaration

end # end module declaration

