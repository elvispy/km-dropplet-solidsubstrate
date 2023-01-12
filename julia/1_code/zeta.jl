using LegendrePolynomials: collectPl
include("./problemConditionStruct.jl");
function zeta(amplitudes::Union{ProblemConditions, Vector{Float64}}; order = length(amplitudes))::Function
    if typeof(amplitudes) <: ProblemConditions; amplitudes = amplitudes.deformation_amplitudes; end
    return (θ::Float64) -> sum(amplitudes .* (collectPl(cos(θ), lmax = order))[1:end])
end