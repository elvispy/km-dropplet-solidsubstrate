function f = zeta_generator(amplitudes, varargin)
    if isstruct(amplitudes); amplitudes = amplitudes.deformation_amplitudes; end
    
    order = length(amplitudes);
    
    if nargin == 1
        warning("I'm deprecated!");
        f = @(theta) arrayfun(@(ang) sum(amplitudes .* legendreP(1:order, cos(ang))), theta, 'UniformOutput', false);
    else
        LP = varargin{1}.LEGENDRE_POLYNOMIALS;
        f = @(theta) arrayfun(@(ang) sum(dot(amplitudes, arrayfun(@(idx) LP{idx}(cos(ang)), 1:order))), theta);
    end
end