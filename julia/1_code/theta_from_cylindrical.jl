using LegendrePolynomials

# LegendrePolys(x::Float64; lmax::Integer) = collectPl(x, lmax = lmax).parent
""" 
    theta_from_cylindrical(r, amplitudes, guess=pi-0.1)::Float64
 Calculates the theta that corresponds to a radius by Newton-Raphson
 Where tfc(theta) = sin(theta) * (R + zeta(theta))

 # Examples 
 ```jldoctest
 julia> theta_from_cylindrical(0, [0, 0])
pi/2
 ```

 """
function theta_from_cylindrical(r, amplitudes; guess = pi-0.1)::Float64
    
    order = length(amplitudes);
    ζ(θ::Float64)   = @. $sum(amplitudes * (collectPl(cos(θ), lmax = order).parent)[2:end]);
    
    # Derivative of the function
    f_prime(θ)      = @. cos(θ) * (1 + ζ(θ)) - sin(θ)^2 * $sum(amplitudes * (collectdnPl(cos(θ), lmax = order, n = 1).parent)[2:end]);
    
    # Function to be minimized
    f_objective(θ)  = @. sin(θ) * (1 + ζ(θ)) - r;

    θ = guess;
    tol_θ = 1e-7;
    n = 1;

    # Newton Method!
    while abs(f_objective(θ)) >= tol_θ && n < 150
        θ = mod(θ - f_objective(θ)./f_prime(θ) - 1e-4, pi/2) + 1e-4 + pi/2; # If solution is close to pi, θ is unstable with mod function (therefore 1e-4 added)
        n = n + 1;
        if n == 50
            θ = 3.14159;
        elseif n == 100
            θ = rand() * pi/2 + pi/2
        end
    end

    return θ
end
