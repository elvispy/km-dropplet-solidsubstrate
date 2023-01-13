
# Struct definitions
struct ProblemConditions
    nb_harmonics::Int64 # Number of legendre polynomials to be used
    deformation_amplitudes::Array{Float64}
    deformation_velocities::Array{Float64}
    pressure_amplitudes::Array{Float64}
    current_time::Float64
    dt::Float64 # The step made to get to this point.
    center_of_mass::Float64 # Current position of the center of mass
    center_of_mass_velocity::Float64
    number_contact_points::Int64
end