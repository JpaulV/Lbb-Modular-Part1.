function plotSpectralRadius(params, spectral_radii, integrationScheme)
    % Plot Spectral Radius vs dt/T for the chosen parameter sets
    % INPUTS:
    %   params         - Struct containing stability parameters (including labels)
    %   spectral_radii - Matrix of spectral radii results
    %   integrationScheme - Numeric value for the integration scheme (1: WBZ-Alpha, 2: Newmark)

    % Set plot title dynamically based on integration scheme
    if integrationScheme == 1
        plotTitle = 'Spectral Radius vs dt / T (WBZ-\alpha and Newmark Methods)';
    elseif integrationScheme == 2
        plotTitle = 'Spectral Radius vs dt / T (Newmark Methods)';
    else
        error('Unsupported integrationScheme provided. Use 1 for WBZ-Alpha or 2 for Newmark.');
    end

    % Plot spectral radii
    figure;
    colors = lines(length(params.alpha_values)); 
    for j = 1:length(params.alpha_values)
        semilogx(params.dt_T_values, spectral_radii(j, :), ...
                 'Color', colors(j, :), ...
                 'LineWidth', 2, ...
                 'DisplayName', params.labels{j}); 
        hold on;
    end

    % Add stability limit and additional plot details
    yline(1, 'k--', 'DisplayName', 'Stability Limit'); 
    xlabel('dt / T', 'FontSize', 14);
    ylabel('Spectral Radius', 'FontSize', 14);
    ylim([0, 1.6]); 
    set(gca, 'YTick', 0:0.2:1.6); 
    title(plotTitle, 'FontSize', 16);
    grid on; 

    % Show the legend and increase its font size
    lgd = legend('show');
    set(lgd, 'FontSize', 14);  % Increase legend text size

    % Also increase the tick labels on the axes
    set(gca, 'FontSize', 14);
end
