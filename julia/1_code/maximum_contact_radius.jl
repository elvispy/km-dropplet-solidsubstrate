using LegendrePolynomials: collectPl, collectdnPl
include("./r_from_spherical.jl"); include("./zeta.jl");
include("./problemConditionStruct.jl");

"""
    theta_max_radius(::Vector{Float64})::Float64

Calculates the radius such that r'(θ) = 0, using newton method
Where r = (1 + ζ(θ))
and 
    dr/dθ = cos(θ) (1 + ζ(θ)) -sin(θ)^2 * ∑Al P'l(cos(θ))

"""
function maximum_contact_radius(amplitudes::Union{ProblemConditions, Vector{Float64}}; guess = pi/2 + pi/4)::Float64
    if typeof(amplitudes) <: ProblemConditions; amplitudes = amplitudes.deformation_amplitudes; end
    order = length(amplitudes);
    ζ = zeta(amplitudes); # sum(amplitudes .* (collectPl(cos(θ), lmax = order))[2:end])
    drdθ(θ::Float64)::Float64 = cos(θ) * (1 + ζ(θ)) - 
        sin(θ)^2 * sum(amplitudes .* (collectdnPl(cos(θ), lmax = order, n = 1))[1:end]);

    dr2dθ2(θ::Float64)::Float64 = - sin(θ) * (1 + ζ(θ)) - 2 * cos(θ) * sin(θ) * 
        sum(amplitudes .* (collectdnPl(cos(θ), lmax = order, n = 1))[1:end]) + 
        sin(θ)^3 * sum(amplitudes .* (collectdnPl(cos(θ), lmax = order, n = 2))[1:end])

    θ = guess;
    tol_θ = 1e-7;
    n = 1;
    # Newton Method!
    while abs(drdθ(θ)) >= tol_θ && n < 150
        θ = mod(θ - drdθ(θ)/dr2dθ2(θ) - 1e-4, pi) + 1e-4; # If solution is close to pi, θ is unstable with mod function (therefore 1e-4 added)
        n = n + 1;
        if n == 50
            θ = 3.14159/2;
        elseif n == 100
            θ = rand()/100 + pi/2
        end
        if n== 149
            println("Hey!")
        end
    end

    return r_from_spherical(θ, amplitudes);

end
