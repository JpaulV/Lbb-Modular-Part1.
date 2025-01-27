% main.m
% Modular Implementation for MDOF 
% Supports WBZ-Alpha, Newmark Methods, and Stability Analysis
% Author: Jerry Paul Varghese, Adcom. Mat.N. 453553
% Date: 14/01/2025
% ReferenceA: AN ALPHA MODIFICATION OF NEWMARKS METHOD(W. L. WOOD, M BOSSAK,0. C. ZIENKIEWICZ )
% Conserving energy and momentum in nonlinear dynamics: A simple implicit time integration scheme (K-J BATHE)

clc; clear;

% Add paths to modular components
addpath('setup', 'solvers', 'postprocessing', 'animation', 'stability');

%% User Inputs with Dropdown Menus

% System Type Selection
choices = {'MDOF Spring mass system'};
[selection, ok] = listdlg('PromptString', 'Select System Type:', ...
                          'SelectionMode', 'single', ...
                          'ListString', choices, ...
                          'ListSize', [300, 100]);
if ~ok
    error('Selection cancelled. Exiting program.');
end
systemType = selection; % 1 for MDOF, 2 for Truss

% Integration Scheme Selection
choices = {'WBZ-Alpha', 'Newmark'};
[selection, ok] = listdlg('PromptString', 'Select Time Integration Scheme:', ...
                          'SelectionMode', 'single', ...
                          'ListString', choices, ...
                          'ListSize', [300, 100]);
if ~ok
    error('Selection cancelled. Exiting program.');
end
integrationScheme = selection; % 1 for WBZ-Alpha, 2 for Newmark
WBZ_alpha_flag = (integrationScheme == 1);
Newmark_flag = (integrationScheme == 2);

% Animation Option
choices = {'Yes', 'No'};
[selection, ok] = listdlg('PromptString', 'Would you like to see an animation?', ...
                          'SelectionMode', 'single', ...
                          'ListString', choices, ...
                          'ListSize', [300, 100]);
if ~ok
    error('Selection cancelled. Exiting program.');
end
animate = (selection == 1); % 1 for Yes, 0 for No

% Stability Analysis Option
choices = {'Yes', 'No'};
[selection, ok] = listdlg('PromptString', 'Perform Stability Analysis?(loading time 5 sec)', ...
                          'SelectionMode', 'single', ...
                          'ListString', choices, ...
                          'ListSize', [300, 100]);
if ~ok
    error('Selection cancelled. Exiting program.');
end
stabilityAnalysisFlag = (selection == 1); % 1 for Yes, 0 for No

%% Simulation Parameters
dt = 0.2681;     % Time step
total_time = 10; % Total simulation time
omega_p = 1.2;   % Prescribed frequency
time = 0:dt:total_time;

% Declare alpha, beta, gamma parameters
alpha = 0;               % WBZ-Alpha: Set alpha (-0.3 to 0), Newmark trapezoidal, uses alpha = 0
beta = (1 - alpha)^2 / 4;   % Derived beta
gamma = 0.5 - alpha;        % Derived gamma

%% System Setup
disp('Setting up the system...');
switch systemType
    case 1  % SDOF System
        disp('System Type: MDOF');
        % Define SDOF system parameters
        prob_params.m2 = 1; prob_params.m3 = 1;
        prob_params.k1 = 1e7; prob_params.k2 = 1;
        prob_params.omega_p = omega_p; prob_params.dt = dt;
        prob_params.total_time = total_time; prob_params.time = time;
        prob_params.M = [prob_params.m2, 0; 0, prob_params.m3];
        prob_params.K = [prob_params.k1 + prob_params.k2, -prob_params.k2; -prob_params.k2, prob_params.k2];
    case 2  % Truss System
        disp('System Type: Truss');
        disp('Truss system implementation is not yet defined. Exiting simulation.');
        return;
end

%% Solver Selection
switch integrationScheme
    case 1
        disp('Using WBZ-Alpha Method...');
        [u, udot, uddot, R] = wbzAlphaSolver(prob_params, alpha, beta, gamma);
    case 2
        disp('Using Newmark Method...');
        [u, udot, uddot, R] = newmarkSolver(prob_params, beta, gamma);
end

%% Post-Processing
if ~isempty(u)
    switch integrationScheme
        case 1
            postProcessing(prob_params, u, udot, uddot, R, 'WBZ-Alpha', alpha, beta, gamma);
        case 2
            postProcessing(prob_params, u, udot, uddot, R, 'Newmark', alpha, beta, gamma);
    end
else
    disp('Post-processing skipped due to invalid solver outputs.');
end

%% Animation (for SDOF system only)
if systemType == 1 && animate == 1
    if ~isempty(u)
        disp('Animating SDOF system...');
        total_anim_time = 5; % Set the desired animation duration in seconds
        animateSDOFSystem(prob_params, u, total_anim_time);
    else
        disp('Animation skipped due to invalid solver outputs.');
    end
end

%% Stability Analysis
if stabilityAnalysisFlag == 1
    disp('Performing Stability Analysis...');

    % Define available parameter sets
    availableSets = struct( ...
        'B1', struct('alpha', -0.1, 'beta', 0.3025, 'gamma', 0.6, ...
                     'label', 'B1: \alpha=-0.1, \beta=0.3025, \gamma=0.6'), ...
        'B2', struct('alpha', -0.1, 'beta', 0.5, 'gamma', 0.6, ...
                     'label', 'B2: \alpha=-0.1, \beta=0.5, \gamma=0.6'), ...
        'B3', struct('alpha', 0.1, 'beta', 0.3025, 'gamma', 0.6, ...
                     'label', 'B3: \alpha=0.1, \beta=0.3025, \gamma=0.6'), ...
        'B4', struct('alpha', 0, 'beta', 0.25, 'gamma', 0.5, ...
                     'label', 'B4: Trapezoidal (\alpha=0, \beta=0.25, \gamma=0.5)'), ...
        'NT', struct('alpha', 0, 'beta', 0.25, 'gamma', 0.5, ...
                     'label', 'NT: Trapezoidal (\beta=1/4, \gamma=1/2)'), ...
        'N2', struct('alpha', 0, 'beta', 0.3, 'gamma', 11/12, ...
                     'label', 'N2: \beta=0.3, \gamma=11/12'), ...
        'N3', struct('alpha', 0, 'beta', 0.5, 'gamma', 0.5, ...
                     'label', 'N3: \beta=0.5, \gamma=0.5'));

    % Perform analysis for WBZ-Alpha methods
    if WBZ_alpha_flag
        disp('Including WBZ-Alpha Methods: B1, B2, B3, B4...');
        paramChoice = {'B1', 'B2', 'B3', 'B4'};

        % Initialize stability_params
        stability_params.alpha_values = [];
        stability_params.beta_values = [];
        stability_params.gamma_values = [];
        stability_params.labels = {};
        stability_params.dt_T_values = logspace(-3, 1, 100);  % dt/T from 10^-3 to 10^1
        stability_params.M = prob_params.M;
        stability_params.K = prob_params.K;
        stability_params.T = 2 * pi / prob_params.omega_p;  % Period of vibration

        % Populate stability_params for WBZ-Alpha
        for i = 1:length(paramChoice)
            setName = paramChoice{i};
            stability_params.alpha_values(end+1) = availableSets.(setName).alpha;
            stability_params.beta_values(end+1) = availableSets.(setName).beta;
            stability_params.gamma_values(end+1) = availableSets.(setName).gamma;
            stability_params.labels{end+1} = availableSets.(setName).label;
        end

        % Call stability analysis and plot for WBZ-Alpha
        spectral_radii = stabilityAnalysis(stability_params);
        plotSpectralRadius(stability_params, spectral_radii, 1); % 1 for WBZ-Alpha

    % Perform analysis for Newmark methods
    elseif Newmark_flag
        disp('Including Newmark Methods: NT, N2, N3...');
        paramChoice = {'NT', 'N2', 'N3'};

        % Initialize stability_params
        stability_params.alpha_values = [];
        stability_params.beta_values = [];
        stability_params.gamma_values = [];
        stability_params.labels = {};
        stability_params.dt_T_values = logspace(-3, 1, 100);  % dt/T from 10^-3 to 10^1
        stability_params.M = prob_params.M;
        stability_params.K = prob_params.K;
        stability_params.T = 2 * pi / prob_params.omega_p;  % Period of vibration

        % Populate stability_params for Newmark
        for i = 1:length(paramChoice)
            setName = paramChoice{i};
            stability_params.alpha_values(end+1) = availableSets.(setName).alpha;
            stability_params.beta_values(end+1) = availableSets.(setName).beta;
            stability_params.gamma_values(end+1) = availableSets.(setName).gamma;
            stability_params.labels{end+1} = availableSets.(setName).label;
        end

        % Call stability analysis and plot for Newmark
        spectral_radii = stabilityAnalysis(stability_params);
        plotSpectralRadius(stability_params, spectral_radii, 2); % 2 for Newmark
    else
        error('No valid method selected. Set WBZ_alpha_flag or Newmark_flag to 1.');
    end
end

%% Completion Message
disp('Simulation completed successfully.');
