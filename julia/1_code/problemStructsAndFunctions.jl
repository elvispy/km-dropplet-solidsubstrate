
module problemStructsAndFunctions
export ProblemConditions, getNextStep, theta_from_cylindrical, LP, sinPoly
using Polynomials
    

include("./problemConditionStruct.jl");

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

sinPoly = Polynomial([(mod(n, 2) == 0 ? 0.0 : ((-1.0)^(n÷2)/factorial(n))) for n = 0:15 ])

end

