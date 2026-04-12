% grid_search_trajectory.m
% Grid search for trajectory-generation hyperparameters (experimental data).
%
% Uses the best nonparametric-simulation hyperparameters (threshold, p) from
% data/opt_hyperparams_nonparamSim.mat to build a truncated data matrix, then
% sweeps over (mu, sigma, lambda) to minimise the trajectory-quality cost
% defined in utils/quantify_function.m.
%
% Output: data/opt_hyperparams_trajOpt.mat  (best_params_trajopt, objective_values)
%
% Prerequisite: data/opt_hyperparams_nonparamSim.mat  (bundled, or run
%               grid_search_nonparametric_sim.m to regenerate).
%
% NOTE: This script may take several hours.
%       Precomputed results are already bundled in data/opt_hyperparams_trajOpt.mat.
%
% Author: Iskandar Khemakhem
% Year:   2024

%% Setup
scriptDir = fileparts(mfilename('fullpath'));
dataDir   = fullfile(scriptDir, 'data');
addpath(genpath(scriptDir));

%% System parameters
dt = 0.05;  % sampling time [s]
m  = 1;     % number of inputs
n  = 6;     % estimated system order
L  = 500;   % trajectory length (time steps)
q  = 5;     % stacked sample size: [ddtheta4, theta1, theta2, theta4, dtheta4]

%% Load experimental data
[H_w, ~, ~, ~] = pageMatrixFromData(L);

%% Boundary conditions  (start: theta4 = 3*pi/8 ≈ 67.5 deg, end: theta4 = pi/2 = 90 deg)
alpha   = 10;   % number of boundary time steps at each end
I_start = 1 : q*alpha;
I_end   = (q*(L-alpha)+1) : q*L;
I_given = [I_start, I_end];

u_start = zeros(alpha, 1);
y_start = [zeros(alpha,1), zeros(alpha,1), (3*pi/8)*ones(alpha,1), zeros(alpha,1)];
w_start = reshape([u_start, y_start]', q*alpha, 1);

u_end = zeros(alpha, 1);
y_end = [zeros(alpha,1), zeros(alpha,1), (pi/2)*ones(alpha,1), zeros(alpha,1)];
w_end = reshape([u_end, y_end]', q*alpha, 1);

w_vec_given = [w_start; w_end];

%% Build truncated data matrix using best simulation hyperparameters
simParams = load(fullfile(dataDir, 'opt_hyperparams_nonparamSim.mat'));
threshold = simParams.best_params.threshold;
p         = simParams.best_params.p;

[~, ~, Perm] = qr(H_w, 'vector');
Perm         = sort(Perm(1:p));
H_w_trunc    = H_w(:, Perm);

[U, S, V] = svd(H_w_trunc, 'econ');
S_vec     = diag(S);
r = find(S_vec >= threshold, 1, 'last');
if isempty(r), r = size(U, 2); end
H_w_trunc = U(:,1:r) * S(1:r,1:r) * V(:,1:r)';
disp(['Rank after truncation: ', num2str(r), '  (p = ', num2str(p), ')']);

%% Precompute constraint matrices (shared across all grid-search trials)
H_w_given = H_w_trunc(I_given, :);

wr = [repmat(w_vec_given(1:q),        alpha,   1); ...
      repmat(w_vec_given(end-q+1:end), L-alpha, 1)];

v = repmat([0, 1, 1, 5, 1], 1, L);
Q = diag(v);
W = diag([ones(q,1); ones(length(I_given)-2*q,1); ones(q,1)]);

end_cond = [w_vec_given(end-2*q+1); w_vec_given(end-2*q+4); w_vec_given(end-q); ...
            w_vec_given(end-q+1);   w_vec_given(end-1);     w_vec_given(end)];
end_mat  = [1, zeros(1,2*q-1); zeros(3), eye(3), zeros(3,4); zeros(2,2*q-2), eye(2)];

c = [w_vec_given(1:2*q); end_cond];
B = [eye(2*q), zeros(2*q, length(w_vec_given)-2*q); ...
     zeros(length(c)-2*q, length(w_vec_given)-2*q), end_mat];

theta_limit = deg2rad(4);
vel_limit   = abs(3 / (30/pi * 243 * (-1/100)));

% M extracts [+theta1; -theta1; +theta2; -theta2; +dtheta4; -dtheta4] at
% every time step from the full q*L trajectory vector.
M = zeros(6*L, q*L);
for ii = 0:(L-1)
    M(6*ii+1, q*ii+2) =  1;   % +theta1
    M(6*ii+2, q*ii+2) = -1;   % -theta1
    M(6*ii+3, q*ii+3) =  1;   % +theta2
    M(6*ii+4, q*ii+3) = -1;   % -theta2
    M(6*ii+5, q*ii+5) =  1;   % +dtheta4
    M(6*ii+6, q*ii+5) = -1;   % -dtheta4
end
d = repmat([theta_limit; theta_limit; theta_limit; theta_limit; vel_limit; vel_limit], L, 1);

% H_tt: finite-difference (total-variation) operator on the channels
% [ddtheta4, theta4, dtheta4], padded with boundary values.
pad_start = diag([w_vec_given(1); w_vec_given(4); w_vec_given(5)]);
pad_end   = diag([w_vec_given(end-q+1); w_vec_given(end-q+4); w_vec_given(end-q+5)]);
row_coef  = sort([(1:q:q*L)'; (4:q:q*L)'; (5:q:q*L)']);

H_tt = [pad_start, zeros(3, p+3); ...
        zeros(length(row_coef), 3), H_w_trunc(row_coef,:), zeros(length(row_coef), 3)] ...
     - [zeros(length(row_coef), 3), H_w_trunc(row_coef,:), zeros(length(row_coef), 3); ...
        zeros(3, p+3), pad_end];

coefs.H_tt      = H_tt;
coefs.H_w_given = H_w_given;
coefs.wr        = wr;
coefs.W         = W;
coefs.B         = B;
coefs.c         = c;
coefs.M         = M;
coefs.d         = d;
coefs.Q         = Q;

%% Grid search
diary(fullfile(scriptDir, 'grid_search_traj_log.txt'));

mu_values     = custom_spaced_values(1e-7, 100, 15, 6);
sigma_values  = custom_spaced_values(1e-7, 100, 15, 7);
lambda_values = custom_spaced_values(1e-7,  10, 15, 7);

[muGrid, sigmaGrid, LambdaGrid] = ndgrid(mu_values, sigma_values, lambda_values);
mu_vector     = muGrid(:);
sigma_vector  = sigmaGrid(:);
lambda_vector = LambdaGrid(:);

objective_values = arrayfun( ...
    @(mu_i, sigma_i, lambda_i) run_traj_generation( ...
        struct('mu', mu_i, 'sigma', sigma_i, 'lambda', lambda_i), ...
        H_w_trunc, w_vec_given, dt, m, n, L, q, coefs), ...
    mu_vector, sigma_vector, lambda_vector);

[~, idx] = min(objective_values);
best_params_trajopt = struct('mu',     mu_vector(idx), ...
                              'sigma',  sigma_vector(idx), ...
                              'lambda', lambda_vector(idx));
disp('Best trajectory-generation parameters:'); disp(best_params_trajopt);
diary off;

save(fullfile(dataDir, 'opt_hyperparams_trajOpt.mat'), ...
     'best_params_trajopt', 'objective_values', ...
     'mu_values', 'sigma_values', 'lambda_values');

fprintf('Saved: %s\n', fullfile(dataDir, 'opt_hyperparams_trajOpt.mat'));

%% ---- Local function ---------------------------------------------------
function cost = run_traj_generation(params, H_w, w_vec_given, dt, m, n, L, q, coefs)
    [~, cost] = generate_trajectory(H_w, w_vec_given, dt, m, n, L, q, coefs, ...
                                     params.lambda, params.mu, params.sigma, 0, 0);
end
