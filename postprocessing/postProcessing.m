function postProcessing(prob_params, u, udot, uddot, R, methodName, alpha, beta, gamma)
    % Handles the visualization of results for WBZ-Alpha and Newmark methods.
    % INPUTS:
    %   prob_params - Struct containing time and problem parameters
    %   u, udot, uddot - Displacement, velocity, and acceleration matrices
    %   R - Reaction force vector
    %   methodName - String indicating the method used ('WBZ-Alpha' or 'Newmark')
    %   alpha, beta, gamma - Parameters used in the method

    time = prob_params.time;

    % Create a new figure for each method to avoid overwriting
    figure('Name', [methodName ' Results'], 'NumberTitle', 'off');
    sgtitle([methodName ' Method Results (α = ', num2str(alpha), ', β = ', num2str(beta), ', γ = ', num2str(gamma), ')']);

    % Plot Displacement
    subplot(2, 3, 1);
    plot(time, u(1, :), 'r-', 'DisplayName', 'Node 2');
    hold on;
    plot(time, u(2, :), 'b-', 'DisplayName', 'Node 3');
    xlabel('Time (s)');
    ylabel('Displacement (m)');
    legend('Location', 'best');
    title('Displacement');
    grid on;

    % Plot Velocity
    subplot(2, 3, 2);
    plot(time, udot(1, :), 'r-', 'DisplayName', 'Node 2');
    hold on;
    plot(time, udot(2, :), 'b-', 'DisplayName', 'Node 3');
    xlabel('Time (s)');
    ylabel('Velocity (m/s)');
    legend('Location', 'best');
    title('Velocity');
    grid on;

    % Plot Acceleration
    subplot(2, 3, 3);
    plot(time, uddot(1, :), 'r-', 'DisplayName', 'Node 2');
    hold on;
    plot(time, uddot(2, :), 'b-', 'DisplayName', 'Node 3');
    xlabel('Time (s)');
    ylabel('Acceleration (m/s^2)');
    legend('Location', 'best');
    title('Acceleration');
    grid on;

    % Plot Reaction Force
    subplot(2, 3, [5, 6]);
    plot(time, R, 'g-', 'DisplayName', 'Reaction at Node 1');
    xlabel('Time (s)');
    ylabel('Reaction Force (N)');
    legend('Location', 'best');
    title('Reaction Force');
    grid on;

    % Display a confirmation message in the command window
    fprintf('Post-processed results for %s method with parameters α = %.3f, β = %.3f, γ = %.3f.\n', ...
            methodName, alpha, beta, gamma);
end
