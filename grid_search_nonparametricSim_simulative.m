% grid_search_nonparametricSim_simulative.m
% Hyperparameter tuning for nonparametric simulation using simulated data.
%
% Generates crane trajectories from the model, constructs Hankel/Page matrices,
% precomputes SVD-truncated data matrices for each (p, threshold) pair, then
% evaluates a grid of (threshold, lambda, p) triples against a reference
% model. Returns the best combination and saves it to data/sim_grid_search.mat.
%
% Precomputed best parameters are already stored in data/sim_grid_search.mat.
% Running this script reproduces that file (may take several hours).
%
% Author: Iskandar Khemakhem
% Year: 2024

%% ----------------------------- Setup -----------------------------------
scriptDir = fileparts(mfilename('fullpath'));
dataDir   = fullfile(scriptDir, 'data');

addpath(genpath(scriptDir));

%% ----------------------------- System parameters -----------------------
dt    = 0.05;   % sampling time [s]
m     = 1;      % number of inputs
n     = 6;      % estimated system order
L     = 300;    % trajectory length (time steps)
q     = 5;      % stacked sample size: [ddtheta4, theta1, theta2, theta4, dtheta4]
noise = false;  % add measurement noise to simulated data

%% ----------------------------- Generate simulated data -----------------
numberOfTraj = 15 * q * L;
[H_w, ud, yd, ~] = generatePageMatrix(numberOfTraj, L, q, dt, noise);

figure
subplot(2,3,1); plot(rad2deg(ud));
xlabel('Time steps ($\Delta t = 0.05\,s$)', 'Interpreter', 'latex');
ylabel('Acceleration ($\mathrm{deg/s}^2$)', 'Interpreter', 'latex');
title('Simulated: input $\ddot{\theta}_4$', 'Interpreter', 'latex')

subplot(2,3,2); plot(rad2deg(yd(:,1)));
xlabel('Time steps ($\Delta t = 0.05\,s$)', 'Interpreter', 'latex');
ylabel('$\theta_1$ [deg]', 'Interpreter', 'latex');
title('Simulated: output $\theta_1$', 'Interpreter', 'latex')

subplot(2,3,3); plot(rad2deg(yd(:,2)));
xlabel('Time steps ($\Delta t = 0.05\,s$)', 'Interpreter', 'latex');
ylabel('$\theta_2$ [deg]', 'Interpreter', 'latex');
title('Simulated: output $\theta_2$', 'Interpreter', 'latex')

subplot(2,3,5); plot(rad2deg(yd(:,3)));
xlabel('Time steps ($\Delta t = 0.05\,s$)', 'Interpreter', 'latex');
ylabel('$\theta_4$ [deg]', 'Interpreter', 'latex');
title('Simulated: output $\theta_4$', 'Interpreter', 'latex')

subplot(2,3,6); plot(rad2deg(yd(:,4)));
xlabel('Time steps ($\Delta t = 0.05\,s$)', 'Interpreter', 'latex');
ylabel('$\dot{\theta}_4$ [deg/s]', 'Interpreter', 'latex');
title('Simulated: output $\dot{\theta}_4$', 'Interpreter', 'latex')

% Verify Willems rank condition on the input
disp(['More data than L?: ', num2str(length(ud) > L)]);
H_u_willems = blkhank(ud, L + n);
disp(['Willems condition satisfied? ', num2str(rank(H_u_willems) >= m*L + n)]);

%% ----------------------------- Build test input sequences -------------
init_cond = 20;
input_sequences = cell(10, 1);
for i = 1:10
    amp  = i / 120;
    vel  = generate_random_input(L - init_cond, amp);
    time = (0:L-init_cond-1) * dt;
    input_sequences{i} = calculateAcceleration(vel, time);
end

figure;
for i = 1:10
    subplot(5, 2, i);
    plot((0:length(input_sequences{i})-1)*dt, input_sequences{i});
    title(sprintf('Acceleration profile %d', i));
    xlabel('Time [s]'); ylabel('[rad/s^2]');
    grid on;
end

%% ----------------------------- Precompute SVD store -------------------
threshold_values = logspace(log10(5e-5), log10(0.05), 20);
lambda_values    = logspace(-5, 0, 20);
p_values         = 1500:500:2500;

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
        input_sequences, dt, m, n, L, q, init_cond), ...
    p_idx, threshold_idx, lambda_vector);

[~, idx]    = min(objective_values);
best_params = struct('threshold', threshold_vector(idx), ...
                     'lambda',    lambda_vector(idx),    ...
                     'p',         p_vector(idx));
disp('Best parameters:'); disp(best_params);

save(fullfile(dataDir, 'sim_grid_search.mat'), ...
     'best_params', 'objective_values', ...
     'threshold_values', 'lambda_values', 'p_values', ...
     'threshold_vector', 'lambda_vector', 'p_vector', ...
     'input_sequences', 'init_cond');

%% ----------------------------- Compare with best params ----------------
simParams = load(fullfile(dataDir, 'sim_grid_search.mat'));
threshold = simParams.best_params.threshold;
lambda    = simParams.best_params.lambda;
p         = simParams.best_params.p;

[~, ~, Perm] = qr(H_w, 'vector');
Perm         = sort(Perm(1:p));
H_w_trunc    = H_w(:, Perm);
[U, S, V]    = svd(H_w_trunc, 'econ');
S_vec        = diag(S);
r = find(S_vec >= threshold, 1, 'last');
if isempty(r), r = size(U, 2); end
H_w_trunc = U(:,1:r) * S(1:r,1:r) * V(:,1:r)';

disp(['Rank after truncation: ', num2str(r), ' (p = ', num2str(p), ')']);

init_cond_idx = 1:q*init_cond;
input_idx     = q*init_cond + 1 : q : q*L;
I_given       = [init_cond_idx, input_idx];

velocity  = generate_random_input(L - init_cond, 0.1);
time_vec  = (0:L-init_cond-1) * dt;
input_seq = calculateAcceleration(velocity, time_vec);

w_vec_given = [repmat([0;0;0;0;0], init_cond, 1); input_seq];
w_dd        = simulate_dd_system(H_w_trunc, I_given, w_vec_given, m, n, L, q, lambda, init_cond);
w_model     = simulate_model(w_vec_given, dt, init_cond, q);

figure;
subplot(2,3,1);
plot(rad2deg(w_dd(1:q:end))); hold on; plot(rad2deg(w_model(1:q:end)));
xlabel('Time steps ($\Delta t=0.05\,s$)', 'Interpreter','latex');
ylabel('Acceleration ($\mathrm{deg/s}^2$)', 'Interpreter','latex');
title('Input $\ddot{\theta}_4$', 'Interpreter','latex');
legend('nonparametric','model');

subplot(2,3,2);
plot(rad2deg(w_dd(2:q:end))); hold on; plot(rad2deg(w_model(2:q:end)));
xlabel('Time steps ($\Delta t=0.05\,s$)', 'Interpreter','latex');
ylabel('$\theta_1$ [deg]', 'Interpreter','latex');
legend('nonparametric','model');

subplot(2,3,3);
plot(rad2deg(w_dd(3:q:end))); hold on; plot(rad2deg(w_model(3:q:end)));
xlabel('Time steps ($\Delta t=0.05\,s$)', 'Interpreter','latex');
ylabel('$\theta_2$ [deg]', 'Interpreter','latex');
legend('nonparametric','model');

subplot(2,3,5);
plot(rad2deg(w_dd(4:q:end))); hold on; plot(rad2deg(w_model(4:q:end)));
xlabel('Time steps ($\Delta t=0.05\,s$)', 'Interpreter','latex');
ylabel('$\theta_4$ [deg]', 'Interpreter','latex');
legend('nonparametric','model');

subplot(2,3,6);
plot(rad2deg(w_dd(5:q:end))); hold on; plot(rad2deg(w_model(5:q:end)));
xlabel('Time steps ($\Delta t=0.05\,s$)', 'Interpreter','latex');
ylabel('$\dot{\theta}_4$ [deg/s]', 'Interpreter','latex');
legend('nonparametric','model');
hold off;

%% ----------------------------- Objective function ----------------------
% Compare nonparametric predictions with model output over 10 test inputs.
%
% Inputs:
%   svd_data       - precomputed low-rank data matrix
%   lambda         - regularisation weight
%   input_sequences - cell(10,1) of acceleration profiles
%   dt, m, n, L, q, init_cond - problem dimensions
%
% Output:
%   cost - accumulated prediction error (lower is better)
function cost = objective_function(svd_data, lambda, input_sequences, dt, m, n, L, q, init_cond) %#ok<INUSD>
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
