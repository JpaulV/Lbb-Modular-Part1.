function animateSDOFSystem(prob_params, u, total_anim_time)
    % Animate the SDOF Spring-Mass System
    % INPUTS:
    %   prob_params      - Struct containing system parameters
    %   u                - Displacement matrix
    %   total_anim_time  - Total time for animation (in seconds)

    % Setup the animation figure
    animFig = figure('Name', 'SDOF System Animation', 'NumberTitle', 'off');
    xlim([-2, 10]);
    ylim([-2, 2]);
    grid on;
    xlabel('X Position (m)');
    ylabel('Y Position (m)');
    title('Spring-Mass System Animation');
    hold on;

    % Plot initial state
    mass_1 = plot(0, 0, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
    text(0, 0.2, 'm1', 'HorizontalAlignment', 'center');
    mass_2 = plot(3, 0, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
    text(3, 0.2, 'm2', 'HorizontalAlignment', 'center');
    mass_3 = plot(6, 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    text(6, 0.2, 'm3', 'HorizontalAlignment', 'center');

    % Spring creation function
    create_spring = @(x1, x2, n_coils) linspace(x1, x2, n_coils * 10);

    % Plot springs
    n_coils = 5;
    x_spring_1 = create_spring(0, 3, n_coils);
    y_spring_1 = 0.1 * sin(2 * pi * n_coils * (x_spring_1 / 3));
    spring_1 = plot(x_spring_1, y_spring_1, 'k-', 'LineWidth', 2);
    text(1.5, 0.5, 'k1', 'HorizontalAlignment', 'center');

    x_spring_2 = create_spring(3, 6, n_coils);
    y_spring_2 = 0.1 * sin(2 * pi * n_coils * (x_spring_2 / 3));
    spring_2 = plot(x_spring_2, y_spring_2, 'k-', 'LineWidth', 2);
    text(4.5, 0.5, 'k2', 'HorizontalAlignment', 'center');

    % Calculate playback speed
    simulation_time = prob_params.total_time; % Total simulation time
    playback_factor = simulation_time / total_anim_time;
    anim_dt = prob_params.dt / playback_factor; % Adjusted time step for animation

    % Animation loop
    for i = 1:length(prob_params.time)
        if ~isvalid(animFig)
            disp('Animation interrupted. Figure was closed.');
            return;
        end

        % Update positions of masses
        set(mass_2, 'XData', 3 + u(1, i));
        set(mass_3, 'XData', 6 + u(2, i));

        % Update springs
        x_spring_1 = create_spring(0, 3 + u(1, i), n_coils);
        y_spring_1 = 0.1 * sin(2 * pi * n_coils * (x_spring_1 / (3 + u(1, i))));
        set(spring_1, 'XData', x_spring_1, 'YData', y_spring_1);

        x_spring_2 = create_spring(3 + u(1, i), 6 + u(2, i), n_coils);
        y_spring_2 = 0.1 * sin(2 * pi * n_coils * (x_spring_2 / (3 + (u(2, i) - u(1, i)))));
        set(spring_2, 'XData', x_spring_2, 'YData', y_spring_2);

        drawnow; % Update the animation
        pause(anim_dt); % Control animation speed
    end

    disp('Animation completed.');
end
