% traj_generation_grid_search.m
% Grid search for trajectory-generation hyperparameters on experimental data.
%
% Loads the bundled experimental crane data, constructs the truncated data
% matrix using the best simulation parameters (from data/sim_grid_search.mat),
% and sweeps over (mu, sigma, lambda) to find the combination that minimises
% the trajectory-quality cost defined in utils/quantify_function.m.
%
% Precomputed best parameters are stored in data/trajopt_grid_search_long.mat.
% Running this script reproduces that file (may take several hours).
%
% Author: Iskandar Khemakhem
% Year: 2024

%% ----------------------------- Setup -----------------------------------
scriptDir = fileparts(mfilename('fullpath'));
dataDir   = fullfile(scriptDir, 'data');

addpath(genpath(scriptDir));

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

% Verify rank conditions
disp(['More data than L?: ',             num2str(length(ud) > L)]);
H_u_willems  = blkhank(ud, L + n);
H_u_berberich = blkhank(ud, L + n + 1);
disp(['Willems condition satisfied? ',   num2str(rank(H_u_willems)  >= m*L + n)]);
disp(['Berberich condition satisfied? ', num2str(rank(H_u_berberich) >= m*(L+n+1))]);

%% ----------------------------- Boundary conditions --------------------
alpha   = 10;   % number of boundary time steps
I_start = 1 : q*alpha;
I_end   = (q*(L-alpha)+1) : q*L;
I_given = [I_start, I_end];

u_start = zeros(alpha, 1);
y_start = [zeros(alpha,1), zeros(alpha,1), (pi/4)*ones(alpha,1), zeros(alpha,1)];
w_start = reshape([u_start, y_start]', q*alpha, 1);

u_end = zeros(alpha, 1);
y_end = [zeros(alpha,1), zeros(alpha,1), (pi/2 - pi/12)*ones(alpha,1), zeros(alpha,1)];
w_end = reshape([u_end, y_end]', q*alpha, 1);

w_vec_given = [w_start; w_end];

disp(['Is given data consistent?: ', ...
      num2str(rank([w_vec_given, H_w(I_given,:)]) == rank(H_w(I_given,:)))]);

%% ----------------------------- Build truncated data matrix -------------
simParams = load(fullfile(dataDir, 'sim_grid_search.mat'));
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

%% ----------------------------- Total-variation matrix ------------------
pad_start = diag([w_vec_given(1); w_vec_given(4); w_vec_given(5)]);
pad_end   = diag([w_vec_given(end-q+1); w_vec_given(end-q+4); w_vec_given(end-q+5)]);
row_coef  = sort([[1:q:q*L]'; [4:q:q*L]'; [5:q:q*L]']);

H_tt = [pad_start, zeros(3, p+3); zeros(length(row_coef), 3), H_w_trunc(row_coef,:), zeros(length(row_coef), 3)] ...
     - [zeros(length(row_coef), 3), H_w_trunc(row_coef,:), zeros(length(row_coef), 3); zeros(3, p+3), pad_end];

%% ----------------------------- Constraint matrices ---------------------
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

M = zeros(6*L, size(H_w_given, 2));
for ii = 0:(L-1)
    M(6*ii+1, q*ii+2) =  1;
    M(6*ii+2, q*ii+2) = -1;
    M(6*ii+3, q*ii+3) =  1;
    M(6*ii+4, q*ii+3) = -1;
    M(6*ii+5, q*ii+5) =  1;
    M(6*ii+6, q*ii+5) = -1;
end
d = repmat([theta_limit; theta_limit; theta_limit; theta_limit; vel_limit; vel_limit], L, 1);

coefs.H_tt      = H_tt;
coefs.H_w_given = H_w_given;
coefs.wr        = wr;
coefs.W         = W;
coefs.B         = B;
coefs.c         = c;
coefs.M         = M;
coefs.d         = d;
coefs.Q         = Q;

%% ----------------------------- Grid search ----------------------------
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

[~, idx]          = min(objective_values);
best_params_trajopt = struct('mu',     mu_vector(idx), ...
                              'sigma',  sigma_vector(idx), ...
                              'lambda', lambda_vector(idx));
disp('Best trajectory-generation parameters:'); disp(best_params_trajopt);
diary off;

save(fullfile(dataDir, 'trajopt_grid_search.mat'), ...
     'best_params_trajopt', 'objective_values', ...
     'mu_values', 'sigma_values', 'lambda_values');

%% ----------------------------- Plot best result ------------------------
trajoptParams = load(fullfile(dataDir, 'trajopt_grid_search_long.mat'));
mu     = trajoptParams.best_params_trajopt.mu;
sigma  = trajoptParams.best_params_trajopt.sigma;
lambda = trajoptParams.best_params_trajopt.lambda;

[w_opt, ~] = generate_trajectory(H_w_trunc, w_vec_given, dt, m, n, L, q, coefs, lambda, mu, sigma, 1, 1);

init_cond           = 20;
I_given_simulation  = [(1:q*init_cond)'; (q*init_cond+1:q:length(w_opt))'];
w_given_simulation  = w_opt(I_given_simulation);
w_dd                = simulate_dd_system(H_w_trunc, I_given_simulation, w_given_simulation, m, n, L, q, lambda, init_cond);
w_model             = simulate_model(w_given_simulation, dt, init_cond, q);

figure;
subplot(2,3,1);
plot(rad2deg(w_dd(1:q:end))); hold on;
plot(rad2deg(w_opt(1:q:end))); plot(rad2deg(w_model(1:q:end)));
xlabel('Time steps ($\Delta t=0.05\,s$)', 'Interpreter','latex');
ylabel('Acceleration ($\mathrm{deg/s}^2$)', 'Interpreter','latex');
title('Input $\ddot{\theta}_4$', 'Interpreter','latex');
legend('nonparametric simulation','trajectory','model');

subplot(2,3,2);
plot(rad2deg(w_dd(2:q:end))); hold on;
plot(rad2deg(w_opt(2:q:end))); plot(rad2deg(w_model(2:q:end)));
xlabel('Time steps ($\Delta t=0.05\,s$)', 'Interpreter','latex');
ylabel('$\theta_1$ [deg]', 'Interpreter','latex');
legend('nonparametric simulation','trajectory','model');

subplot(2,3,3);
plot(rad2deg(w_dd(3:q:end))); hold on;
plot(rad2deg(w_opt(3:q:end))); plot(rad2deg(w_model(3:q:end)));
xlabel('Time steps ($\Delta t=0.05\,s$)', 'Interpreter','latex');
ylabel('$\theta_2$ [deg]', 'Interpreter','latex');
legend('nonparametric simulation','trajectory','model');

subplot(2,3,5);
plot(rad2deg(w_dd(4:q:end))); hold on;
plot(rad2deg(w_opt(4:q:end))); plot(rad2deg(w_model(4:q:end)));
xlabel('Time steps ($\Delta t=0.05\,s$)', 'Interpreter','latex');
ylabel('$\theta_4$ [deg]', 'Interpreter','latex');
legend('nonparametric simulation','trajectory','model');

subplot(2,3,6);
plot(rad2deg(w_dd(5:q:end))); hold on;
plot(rad2deg(w_opt(5:q:end))); plot(rad2deg(w_model(5:q:end)));
xlabel('Time steps ($\Delta t=0.05\,s$)', 'Interpreter','latex');
ylabel('$\dot{\theta}_4$ [deg/s]', 'Interpreter','latex');
legend('nonparametric simulation','trajectory','model');
hold off;

%% ----------------------------- Publication figure ----------------------
time = (1:L) * dt;

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaulttextinterpreter',          'latex');
set(groot, 'defaultLegendInterpreter',        'latex');

fig           = figure('Color', [232 232 232]/255);
subplotHeight = 0.16;
subplotWidth  = 0.80;
subplotLeft   = 0.12;
vSpacing      = 0.18;
c             = 0.80;

positions = [subplotLeft, c,               subplotWidth, subplotHeight;
             subplotLeft, c-vSpacing,      subplotWidth, subplotHeight;
             subplotLeft, c-2*vSpacing,    subplotWidth, subplotHeight;
             subplotLeft, c-3*vSpacing,    subplotWidth, subplotHeight;
             subplotLeft, c-4*vSpacing,    subplotWidth, subplotHeight];

color_ref  = [1,      0.333, 0,      1];
color_opt  = [0,      0.318, 0.620,  1];
color_dd   = [0.624,  0.600, 0.596,  1];
color_model = [0.929, 0.694, 0.125,  1];

ylabels = {'$\ddot{\theta}_4$ [rad/s$^2$]', '$\theta_1$ [rad]', '$\theta_2$ [rad]', ...
           '$\theta_4$ [rad]', '$\dot{\theta}_4$ [rad/s]'};
ylims   = {[-0.1, 0.1], [-0.01, 0.01], [-0.01, 0.01], ...
           [wr(4)-0.05, wr(end-1)+0.05], [-0.05, 0.15]};

for k = 1:5
    hk = subplot(5, 1, k);
    xline(time(alpha),       'Color', [0, 0.447, 0.741, 0.7], 'LineStyle', '-.', 'LineWidth', 1.5); hold on;
    xline(time(end-alpha+1), 'Color', [0, 0.447, 0.741, 0.7], 'LineStyle', '-.', 'LineWidth', 1.5);

    p4   = plot(time, wr(k:q:end),      '--', 'LineWidth', 1.5, 'Color', color_ref);
    p2   = plot(time, w_model(k:q:end), '-',  'LineWidth', 1.5, 'Color', color_model);
    p1   = plot(time, w_dd(k:q:end),    '-',  'LineWidth', 1.5, 'Color', color_dd);
    p3   = plot(time, w_opt(k:q:end),   '-',  'LineWidth', 1.5, 'Color', color_opt);
    sc   = scatter(time(1:15:end), w_opt(k:15*q:end), 'o', ...
                   'MarkerFaceColor', color_opt(1:3), 'MarkerEdgeColor', color_opt(1:3), 'SizeData', 20);
    sc.MarkerFaceAlpha = 0.5;
    hold off;

    ylabel(ylabels{k}, 'FontSize', 14, 'Interpreter', 'latex');
    ylim(ylims{k});
    set(gca, 'FontSize', 14);
    if k < 5, set(gca, 'XTick', []); else, xlabel('Time [s]', 'FontSize', 14, 'Interpreter', 'latex'); end
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

%% ----------------------------- Local helper function ------------------
% Evaluate the trajectory-generation objective for one hyperparameter tuple.
function cost = run_traj_generation(params, H_w, w_vec_given, dt, m, n, L, q, coefs)
    [~, cost] = generate_trajectory(H_w, w_vec_given, dt, m, n, L, q, coefs, ...
                                     params.lambda, params.mu, params.sigma, 0, 0);
end
