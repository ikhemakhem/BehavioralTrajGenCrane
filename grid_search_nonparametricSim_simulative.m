% grid_search_hyperparam_tuning.m
% Hyperparameter tuning (grid search) for nonparametric data-driven simulation
%
% This script generates data, constructs Hankel/Page matrices, precomputes
% SVD-truncated versions of the data matrix for different choices of p and
% singular-value thresholds, evaluates an objective over a grid of
% (threshold, lambda, p) and returns the best combination.
%
% Notes:
% - This is written as a script (not a function). Several helper functions
%   are called but not defined here: generatePageMatrix, blkhank,
%   generate_random_input, calculateAcceleration, simulate_dd_system,
%   simulate_model.
% - The script saves results to 'sim_grid_search.mat' and later reloads them
%   to compare the best found parameters.
% - Variables that should be adapted by the user are grouped near the top.
%
% Author: Iskandar Khemakhem
% Year: 2024

%% ----------------------------- User settings ---------------------------
% basic parameters
dt = 0.05;            % time step [s]
m = 1;                % number of inputs (adjust to your system)
n = 6;                % system order (or expected dimension)
L = 300;              % length of each trajectory (time steps)
q = 5;                % stacked sample size; includes ddtheta4 plus 4 measured outputs

% noise flag for data generation
noise = false;

% how many samples (columns) we request from generatePageMatrix
numberOfTraj = 15 * q * L;

% Construct Page/Hankel-like matrix from generated data
[H_w, ud, yd, wd] = generatePageMatrix(numberOfTraj, L, q, dt, noise);
% H_w: data matrix used to synthesize new trajectories
% ud: input trajectories (samples x ?)
% yd: output trajectories
% wd: full stacked signal (input+outputs) trajectories

%% --------------------------- Quick data checks -------------------------
% Visualize some trajectories: input and selected outputs
figure
subplot(2,3,1);
plot(rad2deg(ud));
title('Data trajectory input ddtheta4')
subplot(2,3,2);
plot(rad2deg(yd(:,1)));
title('Data trajectory output theta1')
subplot(2,3,3);
plot(rad2deg(yd(:,2)));
title('Data trajectory output theta2')
subplot(2,3,5);
plot(rad2deg(yd(:,3)));
title('Data trajectory output theta4')
subplot(2,3,6);
plot(rad2deg(yd(:,4)));
title('Data trajectory output dtheta4')

% Basic verification of Willems/Berberich rank conditions
disp(['Do we have more data than L?: ', num2str(length(ud) > L)]);
H_u_willems = blkhank(ud, L + n);
disp(['Is Willems Lemma condition satisfied? ', num2str(rank(H_u_willems) >= m*L + n)]);

%% ----------------------- Generate random inputs for test ----------------
% fixed initial condition length (in time steps)
init_cond = 20;

% prepare a set of candidate input sequences used by the objective
input_sequences = cell(10,1);
for i = 1:10
    amp = i / 120;                         % amplitude scaling per trial
    vel = generate_random_input(L - init_cond, amp); % generate velocity profile
    time = 0:dt:(L - init_cond - 1) * dt; % matching time vector
    input_sequences{i} = calculateAcceleration(vel, time); % convert to acceleration
end

% plot generated acceleration profiles for inspection
figure;
for i = 1:15
    subplot(5, 2, i);
    plot(time, input_sequences{i});
    title(sprintf('Acceleration Profile %d', i));
    xlabel('Time (s)');
    ylabel('Acceleration (m/s^2)');
    grid on;
end

%% --------------------------- Grid search setup -------------------------
% parameter ranges to sweep
threshold_values = logspace(log10(5e-5), log10(0.05), 20);  % SVD truncation threshold
lambda_values = logspace(-5, 0, 20);                        % regularization parameter
p_values = 1500:500:2500;                                   % number of columns to keep after QR

% keep only unique values (useful for mapping later)
threshold_unique = unique(threshold_values);
p_unique = unique(p_values);

% Preallocate store for SVD-truncated H_w for each unique (p, threshold)
svd_store = cell(numel(p_unique), numel(threshold_unique));

% For each p we first pick the most informative p columns using QR pivoting
% and then compute SVD and apply thresholding to obtain a low-rank approximation.
for j = 1:numel(p_unique)
    p = p_unique(j);

    % QR with column pivoting to pick influential columns of H_w
    [~, ~, P] = qr(H_w, 'vector');
    P = sort(P(1:p));              % choose and sort the first p pivot columns
    H_w_trunc = H_w(:, P);         % truncated data matrix using chosen columns

    for i = 1:numel(threshold_unique)
        threshold = threshold_unique(i);

        % SVD of the truncated matrix
        [U, S, V] = svd(H_w_trunc, 'econ');
        S_vec = diag(S);

        % pick last singular value that is >= threshold
        r = find(S_vec >= threshold, 1, 'last');
        if isempty(r)
            r = size(U, 2); % fallback: use full rank from the SVD
        end

        % reconstruct the truncated (low-rank) approximation
        H_w_trunc_lr = U(:, 1:r) * S(1:r, 1:r) * V(:, 1:r)';

        % store for later reuse in the objective evaluations
        svd_store{j, i} = H_w_trunc_lr;
    end
end

% Create full grid of parameter combinations and flatten to vectors for easy
% mapping to the precomputed SVD store
[ThresholdGrid, LambdaGrid, PGrid] = ndgrid(threshold_values, lambda_values, p_values);
threshold_vector = ThresholdGrid(:);
lambda_vector = LambdaGrid(:);
p_vector = PGrid(:);

% Map each flattened threshold and p to the indices of the precomputed store
[~, threshold_idx] = ismember(threshold_vector, threshold_unique);
[~, p_idx] = ismember(p_vector, p_unique);

%% ------------------------- Objective evaluation ------------------------
% Evaluate the objective for each triple (threshold, lambda, p). We call
% objective_function which runs the data-driven simulation multiple times
% and returns an accumulated error between data-driven and model outputs.
objective_values = arrayfun(@(p_idx, t_idx, l) objective_function(svd_store{p_idx, t_idx}, l, ...
    input_sequences, dt, m, n, L, q, init_cond), ...
    p_idx, threshold_idx, lambda_vector);

% find best (minimum) objective and corresponding parameters
[minValue, idx] = min(objective_values);
best_params = struct('threshold', threshold_vector(idx), 'lambda', lambda_vector(idx), 'p', p_vector(idx));

disp('Best Parameters:');
disp(best_params);

% Save relevant variables for later analysis
save('sim_grid_search.mat', 'best_params', 'objective_values', 'threshold_values', 'lambda_values', 'p_values', ...
     'threshold_vector', 'lambda_vector', 'p_vector', 'input_sequences', 'init_cond');

%% ---------------------- Compare simulations with best params --------------
% reload to ensure reproducible comparison
load('sim_grid_search.mat');
threshold = best_params.threshold;
lambda = best_params.lambda;
p = best_params.p;

% Recompute truncated H_w using the chosen p and threshold (same steps as above)
[~, ~, P] = qr(H_w, 'vector');
P = sort(P(1:p));
H_w_trunc = H_w(:, P);
[U, S, V] = svd(H_w_trunc, 'econ');
S_vec = diag(S);

r = find(S_vec >= threshold, 1, 'last');
if isempty(r), r = size(U,2); end
H_w_trunc = U(:, 1:r) * S(1:r, 1:r) * V(:, 1:r)';

% Indices: initial condition block and subsequent inputs
init_cond_idx = 1:q*init_cond;
input_idx = q*init_cond + 1:q:q*L;
I_given = [init_cond_idx, input_idx];

% Example: create a test input and corresponding stacked vector w_vec_given
velocity = generate_random_input(L - init_cond, 0.1);
time = 0:dt:(L - init_cond - 1) * dt;
input_seq = calculateAcceleration(velocity, time);

w_vec_given = [repmat([0;0;0;0;0], init_cond, 1); input_seq];

% Simulate using data-driven (w_dd) and model (w_model) approaches
w_dd = simulate_dd_system(H_w_trunc, I_given, w_vec_given, m, n, L, q, lambda, init_cond);
w_model = simulate_model(w_vec_given, dt, init_cond, q);

if threshold > 0
    disp(['rank of H_w after truncation: ', num2str(p)]);
    disp(['r ', num2str(r)]);
end

% Plot comparison of key outputs between data-driven and model
figure;
subplot(2,3,1);
plot(rad2deg(w_dd(1:q:end))); hold on;
plot(rad2deg(w_model(1:q:end)));
xlabel('Time steps ($\Delta t = 0.05 \, s$)', 'Interpreter', 'latex');
ylabel('Acceleration ($\mathrm{deg/s}^2$)', 'Interpreter', 'latex');
title('Solution input ddtheta4');
legend('data driven', 'model');

subplot(2,3,2);
plot(rad2deg(w_dd(2:q:end))); hold on;
plot(rad2deg(w_model(2:q:end)));
xlabel('Time steps ($\Delta t = 0.05 \, s$)', 'Interpreter', 'latex');
ylabel('Load sway angle in radial direction ($\mathrm{deg}$)', 'Interpreter', 'latex');
title('Solution output theta1');
legend('data driven', 'model');

subplot(2,3,3);
plot(rad2deg(w_dd(3:q:end))); hold on;
plot(rad2deg(w_model(3:q:end)));
xlabel('Time steps ($\Delta t = 0.05 \, s$)', 'Interpreter', 'latex');
ylabel('Load away angle in tangential direction ($\mathrm{deg}$)', 'Interpreter', 'latex');
title('Solution output theta2');
legend('data driven', 'model');

subplot(2,3,5);
plot(rad2deg(w_dd(4:q:end))); hold on;
plot(rad2deg(w_model(4:q:end)));
xlabel('Time steps ($\Delta t = 0.05 \, s$)', 'Interpreter', 'latex');
ylabel('Horizontal boom angle ($\mathrm{deg}$)', 'Interpreter', 'latex');
title('Solution output theta4');
legend('data driven', 'model');

subplot(2,3,6);
plot(rad2deg(w_dd(5:q:end))); hold on;
plot(rad2deg(w_model(5:q:end)));
xlabel('Time steps ($\Delta t = 0.05 \, s$)', 'Interpreter', 'latex');
ylabel('Horizontal boom velocity ($\mathrm{deg/s}$)', 'Interpreter', 'latex');
title('Solution output dtheta4');
legend('data driven', 'model');
hold off;

%% --------------------- Objective function used in grid search ------------
% Evaluate one hyperparameter tuple by comparing data-driven and model outputs.
function cost = objective_function(svd_data, lambda, input_sequences, dt, m, n, L, q, init_cond)
% cost = objective_function(svd_data, lambda, input_sequences,...)
%   Runs multiple simulations comparing the data-driven predictor (based on
%   the precomputed low-rank data matrix svd_data) with a reference model
%   and returns the accumulated error as the objective to minimize.
%
% Inputs:
% - svd_data: low-rank approximation of the data matrix H_w (columns chosen
%             and truncated by SVD)
% - lambda  : regularization parameter passed to simulate_dd_system
% - input_sequences: cell array of candidate input sequences to test
% - dt, m, n, L, q, init_cond: problem parameters (see top of script)
%
% Output:
% - cost: accumulated error over a fixed number of trials (lower is better)

    total_cost = 0;

    % Precompute index sets for initial conditions and inputs inside w-vector
    init_cond_idx = 1:q*init_cond;
    input_idx = q*init_cond + 1:q:q*L;
    I_given = [init_cond_idx, input_idx];

    % Evaluate over 10 different input trials (this replicates the paper's choice)
    for i = 1:10
        input_seq = input_sequences{i};
        input_seq = input_seq(1:L - init_cond);         % ensure correct length
        w_vec_given = [repmat([0; 0; 0; 0; 0], init_cond, 1); input_seq];

        % predict using data-driven DDS and compare with model
        w_dd = simulate_dd_system(svd_data, I_given, w_vec_given, m, n, L, q, lambda, init_cond);
        w_model = simulate_model(w_vec_given, dt, init_cond, q);

        % Euclidean norm of difference added to total cost
        total_cost = total_cost + norm(w_dd - w_model);
    end

    cost = total_cost;
end
