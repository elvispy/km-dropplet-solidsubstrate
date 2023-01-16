function tt = theta_from_cylindrical(r, amplitudes, PROBLEM_CONSTANTS)
    if isstruct(amplitudes); amplitudes = amplitudes.deformation_amplitudes; end
    
    % order = length(amplitudes);
    zeta = zeta_generator(amplitudes, PROBLEM_CONSTANTS);
    collectdnPl = PROBLEM_CONSTANTS.collectdnPl;
    % zeta(theta::Float64)   = @. $sum(amplitudes * (collectPl(cos(theta), lmax = order).parent)[2:end]);
    
    % Derivative of the function
    f_prime = @(theta) cos(theta) * (1 + zeta(theta)) - sin(theta)^2 * sum(dot(amplitudes, collectdnPl(cos(theta))));
    
    % Function to be minimized00
    f_objective = @(theta) sin(theta) * (1 + zeta(theta)) - r;

    theta = pi - 0.1;
    tol_theta = 1e-7;
    n = 1;

    % Newton Method!
    while abs(f_objective(theta)) >= tol_theta && n < 150
        theta = mod(theta - f_objective(theta)/f_prime(theta) - 1e-4, pi/2) + 1e-4 + pi/2; % If solution is close to pi, theta is unstable with mod function (therefore 1e-4 added)
        n = n + 1;
        if n == 50
            theta = 3.14159;
        elseif n == 100
            theta = rand() * pi/2 + pi/2;
        end
    end

    tt = theta;

end