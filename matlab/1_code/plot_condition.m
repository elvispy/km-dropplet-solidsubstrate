function plot_condition(conditions)
    figure;
    hold on;
    cut = 0.75 * pi;
    sample = [linspace(0, cut, 50), linspace(cut, pi, 100)];
    arrX = sin(sample);
    arrY = cos(sample);
    etas = zeta_generator(conditions.deformation_amplitudes);
    height = conditions.center_of_mass;
    
    EtaX = arrayfun(@(angle) sin(angle) * (1+  etas(angle)), sample);
    EtaY = height + arrayfun(@(angle) cos(angle) .* (1+  etas(angle)), sample);
    
    plot( EtaX,EtaY, 'LineWidth',1.5 , 'Color', [.5 .5 .5]);
    plot(-EtaX,EtaY, 'LineWidth',1.5 , 'Color', [.5 .5 .5]);
    if isstruct(conditions)
        zps = zeta_generator(conditions.pressure_amplitudes);
        ps = @(ang) zps(ang) - sum(conditions.pressure_amplitudes);
        quiver(EtaX, EtaY, -ps(sample) .* arrX, -ps(sample) .* arrY);
    end
    xlim([-1.3, 1.3]);
    ylim([-0.1, 2.5]);
    yline(0, 'r', 'LineWidth', 2);
end