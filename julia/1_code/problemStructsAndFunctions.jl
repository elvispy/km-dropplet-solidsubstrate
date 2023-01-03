module problemStructsAndFunctions
export ProblemConditions, getNextStep
    

# Struct definitions
struct ProblemConditions
    nb_harmonics::Int64 # Number of legendre polynomials to be used
    deformation_amplitudes::Array{Float64}
    velocities_amplitudes::Array{Float64}
    pressure_amplitudes::Array{Float64}
    dt ::Float64
    center_of_mass::Float64 # Current position of the center of mass
    center_of_mass_velocity::Float64
end

include("./getNextStep.jl")

end
