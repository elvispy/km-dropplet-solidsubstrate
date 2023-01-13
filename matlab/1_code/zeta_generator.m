function f = zeta_generator(amplitudes)
    if isstruct(amplitudes); amplitudes = amplitudes.deformation_amplitudes; end
    order = length(theta);
    f = @(theta) sum(amplitudes .* legendreP(1:order, cos(theta)));
end