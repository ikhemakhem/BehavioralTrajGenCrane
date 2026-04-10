% real_data_grid_search_simulation.m
% Grid search for nonparametric-simulation hyperparameters on experimental data.
%
% Loads the bundled experimental crane data and the sim_500 validation set,
% sweeps over (threshold, lambda, p) and returns the best combination by
% comparing nonparametric predictions against real measured trajectories.
%
% Precomputed best parameters are stored in data/sim_grid_search.mat.
% Running this script reproduces that file (may take several hours).
%
% Author: Iskandar Khemakhem
% Year: 2024

%% ----------------------------- Setup -----------------------------------
scriptDir  = fileparts(mfilename('fullpath'));
dataDir    = fullfile(scriptDir, 'data');
sim500Dir  = fullfile(dataDir, 'sim_500');

addpath(genpath(scriptDir));   % ensure all subfolders are on the path

%% ----------------------------- System parameters -----------------------
dt = 0.05;  % sampling time [s]
m  = 1;     % number of inputs
n  = 6;     % estimated system order
L  = 500;   % trajectory length (time steps)
q  = 5;     % stacked sample size: [ddtheta4, theta1, theta2, theta4, dtheta4]

%% ----------------------------- Load experimental data ------------------
[H_w, ud, yd, ~] = pageMatrixFromData(L);

figure
subplot(2,3,1); plot(rad2deg(ud));
xlabel('Time steps ($\Delta t = 0.05\,s$)', 'Interpreter', 'latex');
ylabel('Acceleration ($\mathrm{deg/s}^2$)', 'Interpreter', 'latex');
title('Data: input $\ddot{\theta}_4$', 'Interpreter', 'latex')

subplot(2,3,2); plot(rad2deg(yd(:,1)));
xlabel('Time steps ($\Delta t = 0.05\,s$)', 'Interpreter', 'latex');
ylabel('$\theta_1$ [deg]', 'Interpreter', 'latex');
title('Data: output $\theta_1$', 'Interpreter', 'latex')

subplot(2,3,3); plot(rad2deg(yd(:,2)));
xlabel('Time steps ($\Delta t = 0.05\,s$)', 'Interpreter', 'latex');
ylabel('$\theta_2$ [deg]', 'Interpreter', 'latex');
title('Data: output $\theta_2$', 'Interpreter', 'latex')

subplot(2,3,5); plot(rad2deg(yd(:,3)));
xlabel('Time steps ($\Delta t = 0.05\,s$)', 'Interpreter', 'latex');
ylabel('$\theta_4$ [deg]', 'Interpreter', 'latex');
title('Data: output $\theta_4$', 'Interpreter', 'latex')

subplot(2,3,6); plot(rad2deg(yd(:,4)));
xlabel('Time steps ($\Delta t = 0.05\,s$)', 'Interpreter', 'latex');
ylabel('$\dot{\theta}_4$ [deg/s]', 'Interpreter', 'latex');
title('Data: output $\dot{\theta}_4$', 'Interpreter', 'latex')

% Verify Willems rank condition on the input
disp(['More data than L?: ', num2str(length(ud) > L)]);
H_u_willems = blkhank(ud, L + n);
disp(['Willems condition satisfied? ', num2str(rank(H_u_willems) >= m*L + n)]);

%% ----------------------------- Load validation set ---------------------
init_cond = 20;

input_sequences      = cell(15, 1);
measured_trajectories = cell(15, 1);

for i = 1:15
    inputData  = load(fullfile(sim500Dir, sprintf('input_seq_sim_500_%d.mat',   i)));
    resultsData = load(fullfile(sim500Dir, sprintf('results_seq_sim_500_%d.mat', i)));

    input_sequences{i}       = inputData.acceleration(init_cond+1:end);
    measured_trajectories{i} = organize_measured_trajectory(resultsData);
end

% Visualise acceleration profiles
figure;
for i = 1:15
    subplot(5, 3, i);
    plot(input_sequences{i});
    title(sprintf('Acceleration Profile %d', i));
    xlabel('Time steps'); ylabel('Acceleration [rad/s^2]');
    grid on;
end

%% ----------------------------- Precompute SVD store --------------------
threshold_values = custom_spaced_values(5e-6, 0.005, 5, 4);
lambda_values    = custom_spaced_values(1e-8, 1e-3,  5, 10);
p_values         = 2500:1500:8000;

threshold_unique = unique(threshold_values);
p_unique         = unique(p_values);

svd_store = cell(numel(p_unique), numel(threshold_unique));

for j = 1:numel(p_unique)
    p = p_unique(j);
    [~, ~, Perm] = qr(H_w, 'vector');
    Perm         = sort(Perm(1:p));
    H_w_trunc    = H_w(:, Perm);

    for i = 1:numel(threshold_unique)
        threshold = threshold_unique(i);
        [U, S, V] = svd(H_w_trunc, 'econ');
        S_vec     = diag(S);
        r = find(S_vec >= threshold, 1, 'last');
        if isempty(r), r = size(U, 2); end
        svd_store{j, i} = U(:,1:r) * S(1:r,1:r) * V(:,1:r)';
    end
end

%% ----------------------------- Grid search ----------------------------
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

save(fullfile(dataDir, 'sim_grid_search_real_data_500.mat'), ...
     'best_params', 'objective_values', ...
     'threshold_values', 'lambda_values', 'p_values', ...
     'threshold_vector', 'lambda_vector', 'p_vector', ...
     'input_sequences', 'init_cond');

%% ----------------------------- Compare with best params ----------------
% Recompute truncated H_w for the best (p, threshold) combination.
threshold = best_params.threshold;
lambda    = best_params.lambda;
p         = best_params.p;

[~, ~, Perm] = qr(H_w, 'vector');
Perm         = sort(Perm(1:p));
H_w_trunc    = H_w(:, Perm);
[U, S, V]    = svd(H_w_trunc, 'econ');
S_vec        = diag(S);
r = find(S_vec >= threshold, 1, 'last');
if isempty(r), r = size(U, 2); end
H_w_trunc = U(:,1:r) * S(1:r,1:r) * V(:,1:r)';

init_cond_idx = 1:q*init_cond;
input_idx     = q*init_cond + 5 : q : q*L;
I_given       = [init_cond_idx, input_idx];

for i = 1:20
    inputData   = load(fullfile(sim500Dir, sprintf('input_seq_sim_500_%d.mat',   i)));
    resultsData = load(fullfile(sim500Dir, sprintf('results_seq_sim_500_%d.mat', i)));

    input_seq = inputData.velocity(init_cond+1:end);
    w_meas    = organize_measured_trajectory(resultsData);

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

%% ----------------------------- Objective function ----------------------
% Evaluate one hyperparameter tuple by comparing nonparametric predictions
% against measured trajectories from the sim_500 validation set.
%
% Inputs:
%   svd_data             - low-rank data matrix for the given (p, threshold)
%   lambda               - regularisation weight
%   input_sequences      - cell array of 15 acceleration profiles
%   measured_trajectories - cell array of 15 corresponding measured trajectories
%   m, n, L, q, init_cond - problem dimensions
%
% Output:
%   cost - total prediction error (lower is better)
function cost = objective_function(svd_data, lambda, input_sequences, measured_trajectories, m, n, L, q, init_cond) %#ok<INUSD>
    total_cost    = 0;
    init_cond_idx = 1:q*init_cond;
    input_idx     = q*init_cond + 1 : q : q*L;
    I_given       = [init_cond_idx, input_idx];

    for i = 1:15
        input_seq    = input_sequences{i}(1:L - init_cond);
        measured_traj = measured_trajectories{i};

        w_vec_given = [measured_traj(1:q*init_cond); input_seq];
        w_dd        = simulate_dd_system(svd_data, I_given, w_vec_given, m, n, L, q, lambda, init_cond);
        total_cost  = total_cost + sqrt((w_dd - measured_traj)' * (w_dd - measured_traj));
    end
    cost = total_cost;
end
