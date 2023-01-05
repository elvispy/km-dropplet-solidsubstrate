

module problemStructsAndFunctions
export ProblemConditions, getNextStep, theta_from_cylindrical
    

# Struct definitions
struct ProblemConditions
    nb_harmonics::Int64 # Number of legendre polynomials to be used
    deformation_amplitudes::Array{Float64}
    velocities_amplitudes::Array{Float64}
    pressure_amplitudes::Array{Float64}
    current_time::Float64
    dt::Float64
    center_of_mass::Float64 # Current position of the center of mass
    center_of_mass_velocity::Float64
    number_contact_points::Int64
end

include("./getNextStep.jl");
include("./theta_from_cylindrical.jl");

end

