function r = r_from_spherical(theta, amplitudes, PROBLEM_CONSTANTS)
    d = length(amplitudes);
    LP = PROBLEM_CONSTANTS.LEGENDRE_POLYNOMIALS;
    r = arrayfun(@(angle) ...
        sin(angle) * (1 + sum(dot(amplitudes, ...
        arrayfun(@(idx) LP{idx}(cos(angle)), 1:d)))), theta);
end