using Plots
gr()
t = range(0, 4π, length = 100);
r = range(1, 0, length = 100);

x = cos.(t) .* r;
y = sin.(t) .* r;

@gif for i in eachindex(x)
    scatter(x[i], y[i], lims = (-1, 1), label = "")
end
println("Hey!")
readline()