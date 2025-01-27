% Implements the Newmark time integration method for structural dynamics.

function [u, udot, uddot, R] = newmarkSolver(prob_params, beta, gamma)
    % Extract parameters
    M = prob_params.M;
    K = prob_params.K;
    dt = prob_params.dt;
    time = prob_params.time;
    omega_p = prob_params.omega_p;

    % Initial conditions
    nt = length(time);
    u = zeros(2, nt);    % Displacement at nodes 2 and 3
    udot = zeros(2, nt); % Velocity
    uddot = zeros(2, nt); % Acceleration

    % Initial acceleration
    F_ext = [0; 0];
    uddot(:, 1) = M \ (F_ext - K * u(:, 1));

    % Reaction force at Node 1
    R = zeros(1, nt);

    % Time-stepping loop
    for i = 2:nt
        % Prescribed displacement at Node 1
        u1 = sin(omega_p * time(i));

        % Initial guess for displacements at nodes 2 and 3
        u_guess = u(:, i-1);

        % Newton-Raphson iteration
        tol = 1e-6;
        max_iter = 20;
        for iter = 1:max_iter
            % Residual force vector
            R_res = M * ((1 / (beta * dt^2)) * (u_guess - u(:, i-1)) ...
                        - (1 / (beta * dt)) * udot(:, i-1) ...
                        - (1 / (2 * beta) - 1) * uddot(:, i-1)) ...
                   + K * u_guess - [prob_params.k1 * u1; 0];

            % Convergence check
            if norm(R_res) < tol
                break;
            end

            % Tangent stiffness matrix
            K_tangent = K + (1 / (beta * dt^2)) * M;

            % Incremental displacement
            delta_u = -K_tangent \ R_res;

            % Update guess
            u_guess = u_guess + delta_u;
        end

        % Update displacement, velocity, and acceleration
        u(:, i) = u_guess;
        udot(:, i) = (gamma / (beta * dt)) * (u(:, i) - u(:, i-1)) ...
                    - (gamma / beta - 1) * udot(:, i-1) ...
                    - (gamma / (2 * beta) - 1) * dt * uddot(:, i-1);
        uddot(:, i) = (1 / (beta * dt^2)) * (u(:, i) - u(:, i-1)) ...
                    - (1 / (beta * dt)) * udot(:, i-1) ...
                    - (1 / (2 * beta) - 1) * uddot(:, i-1);

        % Calculate reaction force at Node 1
        R(i) = prob_params.k1 * (u1 - u(1, i));
    end
end
