
include("problemConditionStruct.jl");
using LegendrePolynomials: collectPl
"""
    plot_current_conditions(Union{ProblemConditions, Vector{Float64}})

Plots the current conditions of the 
"""
function plot_current_conditions(data::Union{ProblemConditions, Vector{Float64}}; center_of_mass::Float64 = NaN);

    if typeof(data) <: ProblemConditions; 
        if isnan(center_of_mass)
            center_of_mass = data.center_of_mass
        end
        data = data.deformation_amplitudes; 
    end

    


end