function linterp(xs, ys)
    function(x)
        ((x < xs[1]) || (x > xs[end])) && return NaN
        for i in 1:(length(xs) - 1)
            if xs[i] <= x < xs[i+1]
                l = (x-xs[i]) / (xs[i+1] - xs[i])
                return (1-l) * ys[i] + l * ys[i+1]
            end
        end
        ys[end]
    end
end

function back_euler(F, x0, xn, y0, n)
    h = (xn - x0)/n
    xs = zeros(n+1)
    ys = zeros(n+1)
    xs[1] = x0
    ys[1] = y0
    for i in 1:n
        xs[i + 1] = xs[i] + h
        ## solve y[i+1] = y[i] + h * F(y[i+1], x[i+1])
        ys[i + 1] = find_zero(y -> ys[i] + h * F(y, xs[i + 1]) - y, ys[i]+h)
    end
  linterp(xs, ys)
end

F(y, x; C=1) = sqrt(C/y - 1)
x0, xn, y0 = 0, 1.2, 0
cyc = back_euler(F, x0, xn, y0, 50)
plot(x -> 1 - cyc(x), x0, xn)