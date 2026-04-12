% run_nonparametric_simulation.m
% Nonparametric behavioral simulation validated on experimental crane data.
%
% Builds a truncated data matrix from the bundled training data and
% simulates all 20 sequences in data/validation_data/ using the behavioral
% model, plotting each nonparametric prediction against the measured trajectory.
%
% ---- HOW TO CHOOSE HYPERPARAMETERS ------------------------------------
%   USE_PRECOMPUTED_PARAMS = true  (default)
%       Loads (threshold, lambda, p) from data/opt_hyperparams_nonparamSim_real.mat.
%       Run grid_search_nonparametric_real.m first to generate that file, or
%       set to false and supply values manually (see block below).
%
%   USE_PRECOMPUTED_PARAMS = false
%       Uses the manual values defined in the MANUAL PARAMETERS block.
% -----------------------------------------------------------------------
%
% Author: Iskandar Khemakhem
% Year:   2024

%% Setup
scriptDir = fileparts(mfilename('fullpath'));
dataDir   = fullfile(scriptDir, 'data');
sim500Dir = fullfile(dataDir, 'validation_data');
addpath(genpath(scriptDir));

%% ---- PARAMETER SELECTION ---------------------------------------------
USE_PRECOMPUTED_PARAMS = true;   % true  → load from data/opt_hyperparams_nonparamSim_real.mat
                                  % false → use manual values below

% Manual parameters (only used when USE_PRECOMPUTED_PARAMS = false)
threshold_manual = 1e-4;    % SVD truncation threshold
lambda_manual    = 1e-5;    % l1 regularisation weight
p_manual         = 2500;    % number of QR-selected columns
%% ---------------------------------------------------------------------

%% System parameters
dt = 0.05;  % sampling time [s]
m  = 1;     % number of inputs
n  = 6;     % estimated system order
L  = 500;   % trajectory length (time steps)
q  = 5;     % stacked sample size: [ddtheta4, theta1, theta2, theta4, dtheta4]

%% Load experimental training data
[H_w, ud, ~, ~] = pageMatrixFromData(L);

disp(['More data than L?: ', num2str(length(ud) > L)]);
H_u_willems = blkhank(ud, L + n);
disp(['Willems condition satisfied? ', num2str(rank(H_u_willems) >= m*L + n)]);

%% Load or set hyperparameters
if USE_PRECOMPUTED_PARAMS
    paramFile = fullfile(dataDir, 'opt_hyperparams_nonparamSim_real.mat');
    if ~exist(paramFile, 'file')
        error(['Precomputed parameter file not found:\n  %s\n' ...
               'Run grid_search_nonparametric_real.m first, or set\n' ...
               'USE_PRECOMPUTED_PARAMS = false.'], paramFile);
    end
    simParams = load(paramFile);
    threshold = simParams.best_params.threshold;
    lambda    = simParams.best_params.lambda;
    p         = simParams.best_params.p;
    fprintf('Loaded precomputed params: threshold=%.2e  lambda=%.2e  p=%d\n', ...
            threshold, lambda, p);
else
    % Reached when USE_PRECOMPUTED_PARAMS is set to false above.
    threshold = threshold_manual;
    lambda    = lambda_manual;
    p         = p_manual;
    fprintf('Using manual params: threshold=%.2e  lambda=%.2e  p=%d\n', ...
            threshold, lambda, p);
end

%% Build truncated data matrix
init_cond = 20;

[~, ~, Perm] = qr(H_w, 'vector');
Perm         = sort(Perm(1:p));
H_w_trunc    = H_w(:, Perm);

[U, S, V] = svd(H_w_trunc, 'econ');
S_vec     = diag(S);
r = find(S_vec >= threshold, 1, 'last');
if isempty(r), r = size(U, 2); end
H_w_trunc = U(:,1:r) * S(1:r,1:r) * V(:,1:r)';
disp(['Rank after truncation: ', num2str(r), '  (p = ', num2str(p), ')']);

%% Run nonparametric simulation for all 20 validation sequences
% Known entries: initial conditions (first init_cond steps) + measured
% boom velocity (dtheta4, channel 5) for all remaining steps.
init_cond_idx = 1:q*init_cond;
input_idx     = q*init_cond + 5 : q : q*L;
I_given       = [init_cond_idx, input_idx];

for i = 1:20
    inputData   = load(fullfile(sim500Dir, sprintf('input_seq_sim_500_%d.mat',   i)));
    resultsData = load(fullfile(sim500Dir, sprintf('results_seq_sim_500_%d.mat', i)));

    input_seq   = inputData.velocity(init_cond+1:end);
    w_meas      = organize_measured_trajectory(resultsData);

    w_vec_given = [w_meas(1:q*init_cond); input_seq];
    w_dd        = simulate_dd_system(H_w_trunc, I_given, w_vec_given, m, n, L, q, lambda, init_cond);

    figure;

    subplot(2,3,1);
    plot(rad2deg(w_dd(1:q:end))); hold on; plot(rad2deg(w_meas(1:q:end)));
    xlabel('Time steps ($\Delta t=0.05\,s$)', 'Interpreter','latex');
    ylabel('Acceleration ($\mathrm{deg/s}^2$)', 'Interpreter','latex');
    title(sprintf('Trajectory %d: $\\ddot{\\theta}_4$', i), 'Interpreter','latex');
    legend('nonparametric','measured');

    subplot(2,3,2);
    plot(rad2deg(w_dd(2:q:end))); hold on; plot(rad2deg(w_meas(2:q:end)));
    xlabel('Time steps ($\Delta t=0.05\,s$)', 'Interpreter','latex');
    ylabel('$\theta_1$ [deg]', 'Interpreter','latex');
    legend('nonparametric','measured');

    subplot(2,3,3);
    plot(rad2deg(w_dd(3:q:end))); hold on; plot(rad2deg(w_meas(3:q:end)));
    xlabel('Time steps ($\Delta t=0.05\,s$)', 'Interpreter','latex');
    ylabel('$\theta_2$ [deg]', 'Interpreter','latex');
    legend('nonparametric','measured');

    subplot(2,3,5);
    plot(rad2deg(w_dd(4:q:end))); hold on; plot(rad2deg(w_meas(4:q:end)));
    xlabel('Time steps ($\Delta t=0.05\,s$)', 'Interpreter','latex');
    ylabel('$\theta_4$ [deg]', 'Interpreter','latex');
    legend('nonparametric','measured');

    subplot(2,3,6);
    plot(rad2deg(w_dd(5:q:end))); hold on; plot(rad2deg(w_meas(5:q:end)));
    xlabel('Time steps ($\Delta t=0.05\,s$)', 'Interpreter','latex');
    ylabel('$\dot{\theta}_4$ [deg/s]', 'Interpreter','latex');
    legend('nonparametric','measured');
    hold off;
end
