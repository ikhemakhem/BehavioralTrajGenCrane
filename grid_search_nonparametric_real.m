% grid_search_nonparametric_real.m
% Grid search for nonparametric-simulation hyperparameters (experimental data).
%
% Loads the bundled training data and the validation set, then sweeps over
% (threshold, lambda, p) to find the combination that minimises the
% prediction error on 15 held-out experimental trajectories.
%
% Output: data/opt_hyperparams_nonparamSim_real.mat  (best_params, objective_values)
%
% NOTE: This script may take several hours.
%       Run it once to generate the parameter file needed by
%       run_nonparametric_simulation.m and grid_search_trajectory.m.
%
% Author: Iskandar Khemakhem
% Year:   2024

%% Setup
scriptDir = fileparts(mfilename('fullpath'));
dataDir   = fullfile(scriptDir, 'data');
sim500Dir = fullfile(dataDir, 'validation_data');
addpath(genpath(scriptDir));

%% System parameters
dt = 0.05;  % sampling time [s]
m  = 1;     % number of inputs
n  = 6;     % estimated system order
L  = 500;   % trajectory length (time steps)
q  = 5;     % stacked sample size: [ddtheta4, theta1, theta2, theta4, dtheta4]

%% Load experimental training data
[H_w, ud, ~, ~] = pageMatrixFromData(L);

% Verify Willems rank condition
disp(['More data than L?: ', num2str(length(ud) > L)]);
H_u_willems = blkhank(ud, L + n);
disp(['Willems condition satisfied? ', num2str(rank(H_u_willems) >= m*L + n)]);

%% Load validation set (sequences 1–15 used for the grid search)
init_cond = 20;
input_sequences       = cell(15, 1);
measured_trajectories = cell(15, 1);
for i = 1:15
    inputData   = load(fullfile(sim500Dir, sprintf('input_seq_sim_500_%d.mat',   i)));
    resultsData = load(fullfile(sim500Dir, sprintf('results_seq_sim_500_%d.mat', i)));
    input_sequences{i}       = inputData.acceleration(init_cond+1:end);
    measured_trajectories{i} = organize_measured_trajectory(resultsData);
end

%% Precompute SVD-truncated matrices for each (p, threshold) pair
threshold_values = custom_spaced_values(5e-6, 0.005, 5, 4);
lambda_values    = custom_spaced_values(1e-8, 1e-3,  5, 10);
p_values         = 2500:1500:8000;

threshold_unique = unique(threshold_values);
p_unique         = unique(p_values);

svd_store = cell(numel(p_unique), numel(threshold_unique));
for j = 1:numel(p_unique)
    p_j = p_unique(j);
    [~, ~, Perm] = qr(H_w, 'vector');
    Perm         = sort(Perm(1:p_j));
    H_w_p        = H_w(:, Perm);
    for i = 1:numel(threshold_unique)
        [U, S, V] = svd(H_w_p, 'econ');
        S_vec     = diag(S);
        r = find(S_vec >= threshold_unique(i), 1, 'last');
        if isempty(r), r = size(U, 2); end
        svd_store{j, i} = U(:,1:r) * S(1:r,1:r) * V(:,1:r)';
    end
end

%% Grid search
[ThresholdGrid, LambdaGrid, PGrid] = ndgrid(threshold_values, lambda_values, p_values);
threshold_vector = ThresholdGrid(:);
lambda_vector    = LambdaGrid(:);
p_vector         = PGrid(:);

[~, threshold_idx] = ismember(threshold_vector, threshold_unique);
[~, p_idx]         = ismember(p_vector,         p_unique);

objective_values = arrayfun( ...
    @(pi, ti, li) objective_function(svd_store{pi, ti}, li, ...
        input_sequences, measured_trajectories, m, n, L, q, init_cond), ...
    p_idx, threshold_idx, lambda_vector);

[~, idx]    = min(objective_values);
best_params = struct('threshold', threshold_vector(idx), ...
                     'lambda',    lambda_vector(idx),    ...
                     'p',         p_vector(idx));
disp('Best parameters:'); disp(best_params);

save(fullfile(dataDir, 'opt_hyperparams_nonparamSim_real.mat'), ...
     'best_params', 'objective_values', ...
     'threshold_values', 'lambda_values', 'p_values', ...
     'threshold_vector', 'lambda_vector', 'p_vector', ...
     'input_sequences', 'init_cond');

fprintf('Saved: %s\n', fullfile(dataDir, 'opt_hyperparams_nonparamSim_real.mat'));

%% ---- Local function ---------------------------------------------------
% Evaluate one (threshold, lambda, p) tuple: accumulated prediction error
% between nonparametric simulation and 15 measured trajectories.
function cost = objective_function(svd_data, lambda, input_sequences, measured_trajectories, m, n, L, q, init_cond)
    total_cost    = 0;
    init_cond_idx = 1:q*init_cond;
    input_idx     = q*init_cond + 1 : q : q*L;
    I_given       = [init_cond_idx, input_idx];
    for i = 1:15
        input_seq   = input_sequences{i}(1:L - init_cond);
        w_vec_given = [measured_trajectories{i}(1:q*init_cond); input_seq];
        w_dd        = simulate_dd_system(svd_data, I_given, w_vec_given, m, n, L, q, lambda, init_cond);
        total_cost  = total_cost + norm(w_dd - measured_trajectories{i});
    end
    cost = total_cost;
end
