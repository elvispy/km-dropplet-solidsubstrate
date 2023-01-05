using LegendrePolynomials, OffsetArrays

A = collectPl(0.5, lmax = 5)
amplitudes = [0.0, 1, 2, 3, 4, 5];


function broadcast(f::Any, amplitudes::Vector{Float64}, A::T) where T <: OffsetArrays.OffsetArray 
    result = fill(0.0, length(amplitudes));
    for (ii, jj) in zip(eachindex(amplitudes), eachindex(A))
        result[ii] = f(amplitudes[ii], A[jj]); 
    end

    return result
end


#.+(amplitudes::Vector{Floag64}, A::T) where T <: OffsetArrays.OffsetArray = broadcast(+, amplitudes, A)

println(amplitudes .+ A)