function coefs = project_amplitudes(ff, harmonics_qtt, endpoints, ~) % PROBLEM CONSTANTS used to be here 
    flag = true;
    if flag
        f = @(theta, idx) my_legendre(idx, cos(theta)) .* ff(theta) .* sin(theta);
        try
            f([0 1], 3);
        catch
            f = @(theta, idx) my_legendre(idx, cos(theta))' .* ff(theta) .* sin(theta);
        end
            
        coefs = arrayfun(@(idx) ...
            (2 * idx + 1)/2 * integral(@(theta) f(theta, idx), endpoints(1), endpoints(2), ...
            'AbsTol', 1e-5, 'RelTol', 1e-3), 1:harmonics_qtt);
    else
        if endpoints(2) > PROBLEM_CONSTANTS.nodes(1) || endpoints(1) < PROBLEM_CONSTANTS.nodes(end)
            [PROBLEM_CONSTANTS.nodes, PROBLEM_CONSTANTS.weights] = fclencurt(2^19+1, endpoints(1), endpoints(2));
            warning("Integration vector recalculated");
        end
        loc = 1; n = length(PROBLEM_CONSTANTS.nodes);
        while PROBLEM_CONSTANTS.nodes(loc) > endpoints(1) && 2 * loc < n
            loc = 2 * loc;
        end
        if PROBLEM_CONSTANTS.nodes(loc) > endpoints(1)
            locmax = n;
            locmin = loc;
        else
            locmin = loc/2;
            locmax = loc;
        end

        while locmax - locmin > 1
            idxmed = round((locmax + locmin)/2);
            if PROBLEM_CONSTANTS.nodes(idxmed) > endpoints(1)
                locmin = idxmed;
            else
                locmax = idxmed;
            end
        end
        nodes   = PROBLEM_CONSTANTS.nodes(1:locmax);
        weights = PROBLEM_CONSTANTS.weights(1:locmax);
        if locmax < 1000; warning(srpintf("Too few evaluations for numerical integration: %d", locmax)); end
        
%         coefs = zeros(harmonics_qtt, 1)
%         for idx = 1:harmonics_qtt
%             funcs = my_legendre(idx, cos(theta))
%         end
        coefs = arrayfun(@(idx) ...
            (2*idx+1)/2 * dot(feval(@(theta) my_legendre(idx, cos(theta)) ...
            .* ff(theta) .* sin(theta), nodes), weights), 1:harmonics_qtt);
    end
end