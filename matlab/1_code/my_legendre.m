function y = my_legendre(n, x)
    A = legendre(n, x);
    y = A(1, :);
end