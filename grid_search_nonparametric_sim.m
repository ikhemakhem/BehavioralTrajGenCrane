% grid_search_nonparametric_sim.m
% Grid search for nonparametric-simulation hyperparameters (model-generated data).
%
% Generates crane trajectories from the ODE model, constructs a block Hankel
% data matrix, and sweeps over (threshold, lambda, p) to find the combination
% that minimises the prediction error on held-out model trajectories.
%
% Output: data/opt_hyperparams_nonparamSim.mat  (best_params, objective_values, grid vectors)
%
% NOTE: This script may take several hours.
%       Precomputed results are already bundled in data/opt_hyperparams_nonparamSim.mat.
%
% Author: Iskandar Khemakhem
% Year:   2024

%% Setup
scriptDir = fileparts(mfilename('fullpath'));
dataDir   = fullfile(scriptDir, 'data');
addpath(genpath(scriptDir));

%% System parameters
dt    = 0.05;   % sampling time [s]
m     = 1;      % number of inputs
n     = 6;      % estimated system order
L     = 300;    % trajectory length (time steps)
q     = 5;      % stacked sample size: [ddtheta4, theta1, theta2, theta4, dtheta4]
noise = false;  % add measurement noise to simulated data

%% Generate simulated training data
numberOfTraj = 15 * q * L;
[H_w, ud, ~, ~] = generatePageMatrix(numberOfTraj, L, q, dt, noise);

% Verify Willems rank condition
disp(['More data than L?: ', num2str(length(ud) > L)]);
H_u_willems = blkhank(ud, L + n);
disp(['Willems condition satisfied? ', num2str(rank(H_u_willems) >= m*L + n)]);

%% Build held-out test input sequences
init_cond = 20;
input_sequences = cell(10, 1);
for i = 1:10
    amp  = i / 120;
    vel  = generate_random_input(L - init_cond, amp);
    time = (0:L-init_cond-1) * dt;
    input_sequences{i} = calculateAcceleration(vel, time);
end

%% Precompute SVD-truncated matrices for each (p, threshold) pair
threshold_values = logspace(log10(5e-5), log10(0.05), 20);
lambda_values    = logspace(-5, 0, 20);
p_values         = 1500:500:2500;

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
        input_sequences, dt, m, n, L, q, init_cond), ...
    p_idx, threshold_idx, lambda_vector);

[~, idx]    = min(objective_values);
best_params = struct('threshold', threshold_vector(idx), ...
                     'lambda',    lambda_vector(idx),    ...
                     'p',         p_vector(idx));
disp('Best parameters:'); disp(best_params);

save(fullfile(dataDir, 'opt_hyperparams_nonparamSim.mat'), ...
     'best_params', 'objective_values', ...
     'threshold_values', 'lambda_values', 'p_values', ...
     'threshold_vector', 'lambda_vector', 'p_vector', ...
     'input_sequences', 'init_cond');

fprintf('Saved: %s\n', fullfile(dataDir, 'opt_hyperparams_nonparamSim.mat'));

%% ---- Local function ---------------------------------------------------
% Evaluate one (threshold, lambda, p) tuple: accumulated prediction error
% between nonparametric simulation and model output over 10 test inputs.
function cost = objective_function(svd_data, lambda, input_sequences, dt, m, n, L, q, init_cond)
    total_cost    = 0;
    init_cond_idx = 1:q*init_cond;
    input_idx     = q*init_cond + 1 : q : q*L;
    I_given       = [init_cond_idx, input_idx];
    for i = 1:10
        input_seq   = input_sequences{i}(1:L - init_cond);
        w_vec_given = [repmat([0;0;0;0;0], init_cond, 1); input_seq];
        w_dd        = simulate_dd_system(svd_data, I_given, w_vec_given, m, n, L, q, lambda, init_cond);
        w_model     = simulate_model(w_vec_given, dt, init_cond, q);
        total_cost  = total_cost + norm(w_dd - w_model);
    end
    cost = total_cost;
end
