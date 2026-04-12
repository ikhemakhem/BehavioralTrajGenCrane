% run_trajectory_generation.m
% Optimal crane trajectory generation using the behavioral framework.
%
% Builds a truncated data matrix from the bundled training data and
% generates a minimum-sway optimal boom trajectory from theta4 = 67.5 deg
% to 90 deg.  The result is validated via nonparametric simulation and
% plotted alongside the model-based prediction.
%
% ---- HOW TO CHOOSE HYPERPARAMETERS ------------------------------------
%   USE_PRECOMPUTED_SIM_PARAMS = true  (default)
%       Loads (threshold, p) from data/opt_hyperparams_nonparamSim.mat.
%       Run grid_search_nonparametric_sim.m to regenerate that file.
%
%   USE_PRECOMPUTED_TRAJ_PARAMS = true  (default)
%       Loads (mu, sigma, lambda) from data/opt_hyperparams_trajOpt.mat.
%       Run grid_search_trajectory.m to regenerate that file.
%
%   Set either flag to false and edit the corresponding MANUAL PARAMETERS
%   block to supply your own values.
% -----------------------------------------------------------------------
%
% Author: Iskandar Khemakhem
% Year:   2024

%% Setup
scriptDir = fileparts(mfilename('fullpath'));
dataDir   = fullfile(scriptDir, 'data');
addpath(genpath(scriptDir));

%% ---- PARAMETER SELECTION ---------------------------------------------
USE_PRECOMPUTED_SIM_PARAMS  = true;   % true → load (threshold,p) from opt_hyperparams_nonparamSim.mat
                                       % false → use manual values below
USE_PRECOMPUTED_TRAJ_PARAMS = true;   % true → load (mu,sigma,lambda) from opt_hyperparams_trajOpt.mat
                                       % false → use manual values below

% Manual simulation parameters (used when USE_PRECOMPUTED_SIM_PARAMS = false)
threshold_manual = 1e-4;    % SVD truncation threshold
p_manual         = 2500;    % number of QR-selected columns

% Manual trajectory-generation parameters (used when USE_PRECOMPUTED_TRAJ_PARAMS = false)
mu_manual     = 1.0;     % smoothness (total-variation) penalty weight
sigma_manual  = 1.0;     % reference-tracking penalty weight
lambda_manual = 1e-3;    % l1 sparsity penalty weight
%% ---------------------------------------------------------------------

%% System parameters
dt = 0.05;  % sampling time [s]
m  = 1;     % number of inputs
n  = 6;     % estimated system order
L  = 500;   % trajectory length (time steps)
q  = 5;     % stacked sample size: [ddtheta4, theta1, theta2, theta4, dtheta4]

%% Load experimental data and verify rank conditions
[H_w, ud, ~, ~] = pageMatrixFromData(L);

disp(['More data than L?: ',             num2str(length(ud) > L)]);
H_u_willems   = blkhank(ud, L + n);
H_u_berberich = blkhank(ud, L + n + 1);
disp(['Willems condition satisfied? ',   num2str(rank(H_u_willems)  >= m*L + n)]);
disp(['Berberich condition satisfied? ', num2str(rank(H_u_berberich) >= m*(L+n+1))]);

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

disp(['Is given data consistent?: ', ...
      num2str(rank([w_vec_given, H_w(I_given,:)]) == rank(H_w(I_given,:)))]);

%% Load or set simulation hyperparameters (for data-matrix truncation)
if USE_PRECOMPUTED_SIM_PARAMS
    simParamFile = fullfile(dataDir, 'opt_hyperparams_nonparamSim.mat');
    if ~exist(simParamFile, 'file')
        error(['Simulation parameter file not found:\n  %s\n' ...
               'Run grid_search_nonparametric_sim.m first, or set\n' ...
               'USE_PRECOMPUTED_SIM_PARAMS = false.'], simParamFile);
    end
    simParams = load(simParamFile);
    threshold = simParams.best_params.threshold;
    p         = simParams.best_params.p;
else
    % Reached when USE_PRECOMPUTED_SIM_PARAMS is set to false above.
    threshold = threshold_manual;
    p         = p_manual;
end
fprintf('Sim params: threshold=%.2e  p=%d\n', threshold, p);

%% Build truncated data matrix
[~, ~, Perm] = qr(H_w, 'vector');
Perm         = sort(Perm(1:p));
H_w_trunc    = H_w(:, Perm);

[U, S, V] = svd(H_w_trunc, 'econ');
S_vec     = diag(S);
r = find(S_vec >= threshold, 1, 'last');
if isempty(r), r = size(U, 2); end
H_w_trunc = U(:,1:r) * S(1:r,1:r) * V(:,1:r)';
disp(['Rank after truncation: ', num2str(r), '  (p = ', num2str(p), ')']);

%% Precompute constraint matrices
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
vel_limit   = abs(4 / (30/pi * 243 * (-1/100)));

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

%% Load or set trajectory-generation hyperparameters
if USE_PRECOMPUTED_TRAJ_PARAMS
    trajParamFile = fullfile(dataDir, 'opt_hyperparams_trajOpt.mat');
    if ~exist(trajParamFile, 'file')
        error(['Trajectory parameter file not found:\n  %s\n' ...
               'Run grid_search_trajectory.m first (saves opt_hyperparams_trajOpt.mat),\n' ...
               'or set USE_PRECOMPUTED_TRAJ_PARAMS = false.'], trajParamFile);
    end
    trajParams = load(trajParamFile);
    mu     = trajParams.best_params_trajopt.mu;
    sigma  = trajParams.best_params_trajopt.sigma;
    lambda = trajParams.best_params_trajopt.lambda;
else
    % Reached when USE_PRECOMPUTED_TRAJ_PARAMS is set to false above.
    mu     = mu_manual;
    sigma  = sigma_manual;
    lambda = lambda_manual;
end
fprintf('Traj params: mu=%.2e  sigma=%.2e  lambda=%.2e\n', mu, sigma, lambda);

%% Generate optimal trajectory
[w_opt, ~] = generate_trajectory(H_w_trunc, w_vec_given, dt, m, n, L, q, ...
                                  coefs, lambda, mu, sigma, true, true);

%% Validate via nonparametric simulation and model
init_cond         = 20;
I_given_sim       = [(1:q*init_cond)'; (q*init_cond+1:q:q*L)'];
w_given_sim       = w_opt(I_given_sim);
w_dd              = simulate_dd_system(H_w_trunc, I_given_sim, w_given_sim, m, n, L, q, lambda, init_cond);
w_model           = simulate_model(w_given_sim, dt, init_cond, q);

figure;
ch_labels = {'$\ddot{\theta}_4$ (deg/s$^2$)', '$\theta_1$ (deg)', ...
             '$\theta_2$ (deg)', '$\theta_4$ (deg)', '$\dot{\theta}_4$ (deg/s)'};
sp_idx = [1, 2, 3, 5, 6];
for k = 1:5
    subplot(2, 3, sp_idx(k));
    plot(rad2deg(w_dd(k:q:end)));    hold on;
    plot(rad2deg(w_opt(k:q:end)));
    plot(rad2deg(w_model(k:q:end)));
    xlabel('Time steps ($\Delta t=0.05\,s$)', 'Interpreter','latex');
    ylabel(ch_labels{k}, 'Interpreter','latex');
    legend('nonparametric simulation','trajectory','model');
    hold off;
end

save(fullfile(dataDir, 'opt_traj_result.mat'), 'w_opt', 'w_model', 'w_dd');
fprintf('Saved: %s\n', fullfile(dataDir, 'opt_traj_result.mat'));

%% Results figure
time = (1:L) * dt;

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaulttextinterpreter',          'latex');
set(groot, 'defaultLegendInterpreter',        'latex');

fig           = figure('Color', [232 232 232]/255);
subplotHeight = 0.16;
subplotWidth  = 0.80;
subplotLeft   = 0.12;
vSpacing      = 0.18;
c_pos         = 0.80;

positions = [subplotLeft, c_pos,             subplotWidth, subplotHeight;
             subplotLeft, c_pos-vSpacing,    subplotWidth, subplotHeight;
             subplotLeft, c_pos-2*vSpacing,  subplotWidth, subplotHeight;
             subplotLeft, c_pos-3*vSpacing,  subplotWidth, subplotHeight;
             subplotLeft, c_pos-4*vSpacing,  subplotWidth, subplotHeight];

color_ref   = [1,     0.333, 0,     1];
color_opt   = [0,     0.318, 0.620, 1];
color_dd    = [0.624, 0.600, 0.596, 1];
color_model = [0.929, 0.694, 0.125, 1];

ylabels = {'$\ddot{\theta}_4$ [rad/s$^2$]', '$\theta_1$ [rad]', '$\theta_2$ [rad]', ...
           '$\theta_4$ [rad]', '$\dot{\theta}_4$ [rad/s]'};
ylims   = {[-0.1, 0.1], [-0.01, 0.01], [-0.01, 0.01], ...
           [wr(4)-0.05, wr(end-1)+0.05], [-0.05, 0.15]};

for k = 1:5
    hk = subplot(5, 1, k);
    xline(time(alpha),       'Color', [0, 0.447, 0.741, 0.7], 'LineStyle', '-.', 'LineWidth', 1.5); hold on;
    xline(time(end-alpha+1), 'Color', [0, 0.447, 0.741, 0.7], 'LineStyle', '-.', 'LineWidth', 1.5);

    p4  = plot(time, wr(k:q:end),       '--', 'LineWidth', 1.5, 'Color', color_ref);
    p2  = plot(time, w_model(k:q:end),  '-',  'LineWidth', 1.5, 'Color', color_model);
    p1  = plot(time, w_dd(k:q:end),     '-',  'LineWidth', 1.5, 'Color', color_dd);
    p3  = plot(time, w_opt(k:q:end),    '-',  'LineWidth', 1.5, 'Color', color_opt); 
    sc  = scatter(time(1:15:end), w_opt(k:15*q:end), 'o', ...
                  'MarkerFaceColor', color_opt(1:3), 'MarkerEdgeColor', color_opt(1:3), 'SizeData', 20);
    sc.MarkerFaceAlpha = 0.5;
    hold off;

    ylabel(ylabels{k}, 'FontSize', 14, 'Interpreter', 'latex');
    ylim(ylims{k});
    set(gca, 'FontSize', 14);
    if k < 5
        set(gca, 'XTick', []);
    else
        xlabel('Time [s]', 'FontSize', 14, 'Interpreter', 'latex');
    end
    set(hk, 'Position', positions(k,:));
    box on;

    if k == 1
        dummyLine = plot(nan, nan, 'o-', 'LineWidth', 0.8, 'Color', color_opt);
        legend([p4, dummyLine, p1, p2], ...
               {'Reference', 'Optimal trajectory', 'Nonparametric prediction', 'Model-based prediction'}, ...
               'Location', 'northoutside', 'Orientation', 'horizontal', ...
               'Interpreter', 'latex', 'FontSize', 14, 'Box', 'off');
    end
end

set(fig, 'Units', 'centimeters');
set(fig, 'Position', [0, 0, 22.5, 27]);
set(fig, 'PaperSize', [22.5, 27]);
