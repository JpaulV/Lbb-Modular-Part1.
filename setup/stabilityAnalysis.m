function spectral_radii = stabilityAnalysis(params)
    % Perform stability analysis through spectral radii calculation.
    % INPUT:
    %   params - Struct with fields:
    %       alpha_values  : Vector of alpha values
    %       beta_values   : Vector of beta values
    %       gamma_values  : Vector of gamma values
    %       dt_T_values   : Vector of dt/T ratios
    %       M             : Mass matrix
    %       K             : Stiffness matrix
    %       T             : Period of vibration
    % OUTPUT:
    %   spectral_radii    : Matrix of spectral radii (rows: params, cols: dt/T)
    
    % Initialize spectral radii matrix
    n_params = length(params.alpha_values);
    n_dt_T = length(params.dt_T_values);
    spectral_radii = zeros(n_params, n_dt_T);

    % Loop over dt/T values
    for i = 1:n_dt_T
        dt_T = params.dt_T_values(i);
        dt = dt_T * params.T;  % Time step based on dt/T
        
        % Loop over parameter combinations
        for j = 1:n_params
            alpha = params.alpha_values(j);
            beta = params.beta_values(j);
            gamma = params.gamma_values(j);

            % Construct matrices H1 and H0
            H1 = params.M + (beta + alpha) * dt^2 * params.K;
            H0 = params.M - ((0.5 - beta) + alpha) * dt^2 * params.K;

            % Check if H1 is invertible
            if det(H1) == 0
                error('Matrix H1 is singular and cannot be inverted at dt/T = %f', dt_T);
            end

            % Compute amplification matrix A
            A = H1 \ H0;

            % Calculate eigenvalues of A
            eigenvalues = eig(A);

            % Calculate spectral radius
            spectral_radii(j, i) = max(abs(eigenvalues));
        end
    end
end
