
module problemStructsAndFunctions
export ProblemConditions, getNextStep, theta_from_cylindrical, LP, sinPoly, integrate_poly
using Polynomials
    

# Struct definitions
struct ProblemConditions
    nb_harmonics::Int64 # Number of legendre polynomials to be used
    deformation_amplitudes::Array{Float64}
    velocities_amplitudes::Array{Float64}
    pressure_amplitudes::Array{Float64}
    current_time::Float64
    dt::Float64 # The step made to get to this point.
    center_of_mass::Float64 # Current position of the center of mass
    center_of_mass_velocity::Float64
    number_contact_points::Int64
end

include("./getNextStep.jl");
include("./theta_from_cylindrical.jl");


"""
    LP(::Int)::Vector{Polynomial{Rational{BigInt}, :x}}

Returns a list of the first n Legendre Polynomials (from o1 to n)
"""
function LP(n::Int)::Vector{Polynomial{Rational{BigInt}, :x}}
    res = fill(Polynomial([Rational(BigInt(1))]), (n+1, ));
    x = Polynomial([Rational{BigInt}(0), Rational{BigInt}(1)])
    for ii = 1:n
        if ii == 1
            res[ii+1] = Polynomial([0.0, 1.0]);
        else
            res[ii+1] = ((2*ii-1) * x * res[ii] - (ii-1) * res[ii-1])/ii
        end
    end
    return res[2:end]
end

# This function calculates the primitive of P/(x^-d):  (P::Polynomial) -> ∫P(x) * x^(-d) dx
function integrate_poly(p::Polynomial{T, :x}, d::Int = 1) where T<:Real 
    y = integrate(Polynomial([p.coeffs[ii] for ii = (1+d):(degree(p)+1)]));
    c = p.coeffs[1:d];
    f(x) = sum(c[1:(end-1)] .* [-1/(m*x^m) for m = (d-1):-1:1]) + c[end] * log(abs(x))
    return (t) -> y(t) + f(t)
end

sinPoly = Polynomial([(mod(n, 2) == 0 ? 0.0 : ((-1.0)^(n÷2)/factorial(n))) for n = 0:15 ])
end

#using .problemStructsAndFunctions
#using Polynomials

#println(integrate_poly(Polynomial([0.0, 1.0]), 1)(-1.91287))


