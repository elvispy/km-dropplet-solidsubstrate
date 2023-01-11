using LegendrePolynomials

function r_from_spherical(θ::Float64, amplitudes::Vector{Float64})::Float64
    return sin(θ) * (1 + sum(amplitudes .* (collectPl(cos(θ), lmax = length(amplitudes)).parent)[2:end]));
end