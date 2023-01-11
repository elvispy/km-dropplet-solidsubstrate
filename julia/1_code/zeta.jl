using LegendrePolynomials: collectPl
function zeta(amplitudes; order = length(amplitudes))::Function
    return (θ::Float64) -> sum(amplitudes .* (collectPl(cos(θ), lmax = order))[1:end])
end