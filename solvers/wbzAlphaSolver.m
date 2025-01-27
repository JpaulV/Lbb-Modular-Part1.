% wbzAlphaSolver.m
% Implements the WBZ-Alpha time integration scheme.

function [u, udot, uddot, R] = wbzAlphaSolver(prob_params, alpha, beta, gamma)
    % Extract parameters
    M = prob_params.M;
    K = prob_params.K;
    dt = prob_params.dt;
    time = prob_params.time;
    omega_p = prob_params.omega_p;

    % Initial conditions
    nt = length(time);
    u = zeros(2, nt);
    udot = zeros(2, nt);
    uddot = zeros(2, nt);

    % Initial acceleration
    F_ext = [0; 0];
    uddot(:, 1) = M \ (F_ext - K * u(:, 1));

    % Reaction force at Node 1
    R = zeros(1, nt);

    % Time-stepping loop
    for i = 2:nt
        % Prescribed displacement at Node 1
        u1 = sin(omega_p * time(i));

        % Effective external force
        F_eff = [prob_params.k1 * u1; 0] * (1 + alpha);

        % Predictor step
        u(:, i) = u(:, i-1) + dt * udot(:, i-1) + (0.5 - beta) * dt^2 * uddot(:, i-1);
        udot(:, i) = udot(:, i-1) + (1 - gamma) * dt * uddot(:, i-1);

        % Residual
        R_res = M * uddot(:, i) + K * u(:, i) - F_eff;

        % Newton-Raphson Iteration
        tol = 1e-6;
        max_iter = 20;
        for iter = 1:max_iter
            K_tangent = K + (1 / (beta * dt^2)) * M;
            delta_u = -K_tangent \ R_res;
            u(:, i) = u(:, i) + delta_u;
            R_res = M * ((1 / (beta * dt^2)) * (u(:, i) - u(:, i-1)) ...
                         - (1 / (beta * dt)) * udot(:, i-1) ...
                         - (1 / (2 * beta) - 1) * uddot(:, i-1)) ...
                    + K * u(:, i) - F_eff;
            if norm(R_res) < tol
                break;
            end
        end

        % Update velocity and acceleration
        uddot(:, i) = (1 / (beta * dt^2)) * (u(:, i) - u(:, i-1)) ...
                    - (1 / (beta * dt)) * udot(:, i-1) ...
                    - (1 / (2 * beta) - 1) * uddot(:, i-1);
        udot(:, i) = udot(:, i-1) + dt * ((1 - gamma) * uddot(:, i-1) + gamma * uddot(:, i));

        % Reaction force at Node 1
        R(i) = prob_params.k1 * (u1 - u(1, i));
    end
end
