module problemStructsAndFunctions
export ProblemConditions
    

# Struct definitions
struct ProblemConditions
    nb_harmonics::Int64
    amplitudes::Array{Float64}
    amplitudes_velocities::Array{Float64}
    pressure_disribution::Array{Float64}
    dt ::Float64
end

end
