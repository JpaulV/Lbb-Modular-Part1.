% setupSystem.m
% Sets up the system parameters and returns structures.

function [prob_params, matProps] = setupSystem(systemType, omega_p, dt, total_time, time)
    switch systemType
        case 1
            % Example: MDOF system
            prob_params.m2 = 1;   % Mass at Node 2
            prob_params.m3 = 1;   % Mass at Node 3
            prob_params.k1 = 1e7; % Stiffness between Node 1 and Node 2
            prob_params.k2 = 1;   % Stiffness between Node 2 and Node 3
            prob_params.omega_p = omega_p; % Prescribed frequency
            prob_params.dt = dt;  % Time step size
            prob_params.total_time = total_time; % Total simulation time
            prob_params.time = time; % Time vector

            % Mass and stiffness matrices
            prob_params.M = [prob_params.m2, 0; 
                              0, prob_params.m3];
            prob_params.K = [prob_params.k1 + prob_params.k2, -prob_params.k2; 
                            -prob_params.k2, prob_params.k2];

        otherwise
            error('Unsupported system type.');
    end

    % Material properties (optional)
    matProps = struct();  % Empty for now
end
