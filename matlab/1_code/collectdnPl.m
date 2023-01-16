function res = collectdnPl(x, nmax, varargin) % TODO: Change to cell of function handles created on compile time. It is more efficient
    res = 1;    
    error("Hey, I'm deprecated!");
    if nargin == 2 || varargin{1} == 1
        %lp = x * ones(nmax, 1);
        res = ones(1, nmax);
        %val   = x;
        %valdn = 1;
        for ii = 2:nmax
            res(ii) = ii * legendreP(ii-1, x) + x * res(ii-1);
        end
    elseif varargin{1} == 2
        v = collectdnPl(x, nmax, 1);
        
        res = zeros(1, nmax);
        for ii = 2:nmax
           res(ii) = (ii+1) * v(ii-1) + x * res(ii-1);
        end
    end
end