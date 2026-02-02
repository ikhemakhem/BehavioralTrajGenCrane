%%
%load("sim_grid_search.mat");
%%
% Datd riven trajectory Generation
dt = 0.05; % time step
m = 1;
n = 6;
L = 500;
q = 5; % input + output = size(w(t))

%% 1) Generate data
%buffer = 1000;  % 1193 =  max buffer -1
%noise = true;
%numberOfTraj = 15*q*L;
%[H_w, ud, yd, wd] = generatePageMatrix(numberOfTraj, L, q, dt, noise);
[H_w, ud, yd, wd] = pageMatrixData(L);

figure
subplot(2,3,1);
plot(rad2deg(ud));
xlabel('Time steps ($\Delta t = 0.05 \, s$)', 'Interpreter', 'latex');
ylabel('Acceleration ($\mathrm{deg/s}^2$)', 'Interpreter', 'latex');
title('Data trajectory input ddtheta4')
subplot(2,3,2);
plot(rad2deg(yd(:,1)));
xlabel('Time steps ($\Delta t = 0.05 \, s$)', 'Interpreter', 'latex');
ylabel('Load sway angle in radial direction ($\mathrm{deg}$)', 'Interpreter', 'latex');
title('Data trajectory output theta1')
subplot(2,3,3);
plot(rad2deg(yd(:,2)));
xlabel('Time steps ($\Delta t = 0.05 \, s$)', 'Interpreter', 'latex');
ylabel('Load away angle in tangential direction ($\mathrm{deg}$)', 'Interpreter', 'latex');
title('Data trajectory output theta2')
subplot(2,3,5);
plot(rad2deg(yd(:,3)));
xlabel('Time steps ($\Delta t = 0.05 \, s$)', 'Interpreter', 'latex');
ylabel('Horizontal boom angle ($\mathrm{deg}$)', 'Interpreter', 'latex');
title('Data trajectory output theta4')
subplot(2,3,6);
plot(rad2deg(yd(:,4)));
xlabel('Time steps ($\Delta t = 0.05 \, s$)', 'Interpreter', 'latex');
ylabel('Horizontal boom velocity ($\mathrm{deg/s}$)', 'Interpreter', 'latex');
title('Data trajectory output dtheta4')

% Verify
disp(['Do we have more data than L?: ', num2str(length(ud)> L)]);
% Test input condition
H_u_willems = blkhank(ud, L+n);
disp(['Is Willems Lemma condition satisfied? ', num2str(rank(H_u_willems) >= m*L+n)]);
H_u_berberich = blkhank(ud, L+n+1);
disp(['Is Berberich condition satisfied? ', num2str(rank(H_u_berberich) >= m*(L+n+1))]);

%% 2) Define given data
% init and end condition length
alpha = 10;
I_start = 1:q*alpha;
I_end = q*(L-alpha)+1 : q*L;
I_given = [I_start, I_end];

u_start = zeros(alpha,1);
y_start = [zeros(alpha,1), zeros(alpha,1), pi/4*ones(alpha,1), zeros(alpha,1)];
w_start = [u_start, y_start];
w_start = reshape(w_start', 1, q*alpha)';

u_end = zeros(alpha,1);
y_end = [zeros(alpha,1), zeros(alpha,1), (pi/2-pi/12)*ones(alpha,1),zeros(alpha,1)];
w_end = [u_end, y_end];
w_end = reshape(w_end', 1, q*alpha)';

H_w_given = H_w(I_given,:);

w_vec_given = [w_start; w_end];

% Test consistency
disp(['Is given data consistent?: ', num2str(rank([w_vec_given, H_w_given]) == rank(H_w_given))]);

%% Load the results from simulation grid search
load('sim_grid_search.mat');
threshold = best_params.threshold; %best_params.threshold;
%lambda = best_params.lambda;
p = best_params.p;

% Perform QR decomposition on H_w
[~,~,P] = qr(H_w, 'vector');
P = sort(P(1:p));  % Select the first p columns in sorted order
H_w_trunc = H_w(:, P);

[U, S, V] = svd(H_w_trunc);
S_vec = diag(S);
r = find(S_vec >= threshold, 1, 'last'); % Last significant singular value
if isempty(r)
    r = size(U, 2); % Use all vectors if threshold is too small
end
H_w_trunc = U(:, 1:r)*S(1:r, 1:r)*V(:, 1:r)';

%% Total Variation matrix
pad_start = diag([w_vec_given(1); w_vec_given(4); w_vec_given(5)]);
pad_end = diag([w_vec_given(end-q+1); w_vec_given(end-q+4); w_vec_given(end-q+5)]);
row_coef = [[1:q:q*L]'; [4:q:q*L]'; [5:q:q*L]'];
row_coef = sort(row_coef);

H_tt = [pad_start, zeros(3,p+3); zeros(length(row_coef),3), H_w_trunc(row_coef,:),zeros(length(row_coef),3)]...
    - [zeros(length(row_coef),3), H_w_trunc(row_coef,:), zeros(length(row_coef),3); zeros(3,p+3), pad_end];

%% Constraints
H_w_given = H_w_trunc(I_given,:);
    
wr = [repmat(w_vec_given(1:q), alpha, 1); repmat(w_vec_given(end-q+1:end), L-alpha, 1)];
%gr = H_w_trunc\wr;

v = [0, 1, 1, 5, 1];
v = repmat(v,1,L);
Q = diag(v);

W = diag([1*ones(q,1);ones(length(I_given)-2*q,1);1*ones(q,1)]);

end_cond = [w_vec_given(end-2*q+1); w_vec_given(end-2*q+4); w_vec_given(end-q); w_vec_given(end-q+1); w_vec_given(end-1); w_vec_given(end)];
end_mat = [1, zeros(1,2*q-1); zeros(3), eye(3), zeros(3,4); zeros(2,2*q-2), eye(2)];

c = [w_vec_given(1:2*q); end_cond];
B = [eye(2*q), zeros(2*q,length(w_vec_given)-2*q);zeros(length(c)-2*q,length(w_vec_given)-2*q), end_mat];

theta_limit = deg2rad(4); %[rad]
% put limit on the velocity
vel_limit = abs(3/(30/pi*243*(-1/100)));
% Initialize the selection matrix with zeros
M = zeros(6*L, length(size(H_w_given, 2)));
% Fill the selection matrix
for ii = 0:(L - 1)
    M(6*ii + 1, q*ii + 2) = 1;
    M(6*ii + 2, q*ii + 2) = -1;
    M(6*ii + 3, q*ii + 3) = 1;
    M(6*ii + 4, q*ii + 3) = -1;
    M(6*ii + 5, q*ii + 5) = 1;
    M(6*ii + 6, q*ii + 5) = -1;
end
d = repmat([theta_limit; theta_limit; theta_limit; theta_limit; vel_limit; vel_limit], L,1);

coefs.H_tt = H_tt;
coefs.H_w_given = H_w_given;
coefs.wr = wr;
coefs.W = W;
coefs.B = B;
coefs.c = c;
coefs.M = M;
coefs.d = d;
coefs.Q = Q;

%% Grid Search optimization over hyperparameters
% start recording command window for logging
diary('grid_search_log_real_data.txt'); % Start recording
% Define the hyperparameter space
% Define parameter ranges
mu_values = custom_spaced_values(1e-7, 100, 15, 6);  %15
sigma_values = custom_spaced_values(1e-7, 100, 15, 7);
lambda_values = custom_spaced_values(1e-7, 10, 15, 7);
% mu_values = [1e-7, 1e-4, 5e-2,  1, 10, 25, 100];
% sigma_values = [1e-7,  5e-5, 5e-2,  1, 10, 25, 100];
% lambda_values = [1e-7, 5e-5, 5e-2, 1, 10, 25, 100];

% Generate grid of all combinations
[muGrid, sigmaGrid, LambdaGrid] = ndgrid(mu_values, sigma_values, lambda_values);

% Flatten grids to enable vectorized evaluation
mu_vector = muGrid(:);
sigma_vector = sigmaGrid(:);
lambda_vector = LambdaGrid(:);

% Assuming H_w, I_given, w_vec_given, dt, m, n, L, q, alpha, lambda are defined outside this snippet
objective_values = arrayfun(@(t, s, l) run_traj_generation(struct('mu', t, 'sigma', s, 'lambda', l), ...
    H_w_trunc, w_vec_given, dt, m, n, L, q, coefs), mu_vector, sigma_vector, lambda_vector);

% Find the index of the minimum objective value
[minValue, idx] = min(objective_values);

% Extract the best parameters
best_params_trajopt = struct('mu', mu_vector(idx), 'sigma', sigma_vector(idx), 'lambda', lambda_vector(idx));

% Display the best parameters
disp('Best Parameters:');
disp(best_params_trajopt);
diary off;             % Stop recording
save('trajopt_grid_search_real_data.mat', 'best_params_trajopt', 'objective_values', 'mu_values', 'sigma_values','lambda_values');
%% Plot best result
load('trajopt_grid_search_long.mat');
mu = best_params_trajopt.mu; %1e-4;0
sigma = best_params_trajopt.sigma;
lambda = best_params_trajopt.lambda;
[w_opt, cost] = generate_trajectory(H_w_trunc, w_vec_given, dt, m, n, L, q, coefs, lambda, mu, sigma, 1, 1);

%%
init_cond = 20;
I_given_simulation = [[1:q*init_cond]';[q*init_cond+1:q:length(w_opt)]'];
w_given_simulation = w_opt(I_given_simulation);
w_dd = simulate_dd_system(H_w_trunc, I_given_simulation, w_given_simulation, m, n, L, q, lambda, init_cond);
w_model = simulate_model(w_given_simulation, dt, init_cond, q);

%%
figure;
subplot(2,3,1);
plot(rad2deg(w_dd(1:q:end)));
hold on;
plot(rad2deg(w_opt(1:q:end)));
plot(rad2deg(w_model(1:q:end)));
xlabel('Time steps ($\Delta t = 0.05 \, s$)', 'Interpreter', 'latex');
ylabel('Acceleration ($\mathrm{deg/s}^2$)', 'Interpreter', 'latex');
title('Input ddtheta4');
legend('dd simulation', 'trajectory', 'model');
subplot(2,3,2);
plot(rad2deg(w_dd(2:q:end)));
hold on;
plot(rad2deg(w_opt(2:q:end)));
plot(rad2deg(w_model(2:q:end)));
xlabel('Time steps ($\Delta t = 0.05 \, s$)', 'Interpreter', 'latex');
ylabel('Load sway angle in radial direction ($\mathrm{deg}$)', 'Interpreter', 'latex');
title('Output theta1');
legend('dd simulation', 'trajectory', 'model');
subplot(2,3,3);
plot(rad2deg(w_dd(3:q:end)));
hold on;
plot(rad2deg(w_opt(3:q:end)));
plot(rad2deg(w_model(3:q:end)));
xlabel('Time steps ($\Delta t = 0.05 \, s$)', 'Interpreter', 'latex');
ylabel('Load away angle in tangential direction ($\mathrm{deg}$)', 'Interpreter', 'latex');
title('Output theta2');
legend('dd simulation', 'trajectory', 'model');
subplot(2,3,5);
plot(rad2deg(w_dd(4:q:end)));
hold on;
plot(rad2deg(w_opt(4:q:end)));
plot(rad2deg(w_model(4:q:end)));
xlabel('Time steps ($\Delta t = 0.05 \, s$)', 'Interpreter', 'latex');
ylabel('Horizontal boom angle ($\mathrm{deg}$)', 'Interpreter', 'latex');
title('Output theta4');
legend('dd simulation', 'trajectory', 'model');
hold off;
subplot(2,3,6);
plot(rad2deg(w_dd(5:q:end)));
hold on;
plot(rad2deg(w_opt(5:q:end)));
plot(rad2deg(w_model(5:q:end)));
xlabel('Time steps ($\Delta t = 0.05 \, s$)', 'Interpreter', 'latex');
ylabel('Horizontal boom velocity ($\mathrm{deg/s}$)', 'Interpreter', 'latex');
title('Output dtheta4');
legend('dd simulation', 'trajectory', 'model');
hold off;

%%
% Assuming w_dd and w_model are column vectors and q is a positive integer
    
% Calculate the time vector
time = 0.05:0.05:length(w_dd)*0.2*0.05;

% Set default LaTeX interpreters for labels
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaulttextinterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

% Create a new figure with a specific background color
fig = figure('Color', [232/255, 232/255, 232/255]);
subplotHeight = 0.16; % height of each subplot
subplotWidth = 0.8;   % width of each subplot
subplotLeft = 0.12;    % left position of each subplot
subplotVerticalSpacing = 0.18; % spacing between subplots
c = 0.8;
% Predefine subplot positions
positions = [
    subplotLeft, c, subplotWidth, subplotHeight;
    subplotLeft, c-subplotVerticalSpacing, subplotWidth, subplotHeight;
    subplotLeft, c-2*subplotVerticalSpacing, subplotWidth, subplotHeight;
    subplotLeft, c-3*subplotVerticalSpacing, subplotWidth, subplotHeight;
    subplotLeft, c-4*subplotVerticalSpacing, subplotWidth, subplotHeight;
];

color1 = [1, 0.333, 0, 1]; % Yellow-Orange 
color2 = [0.0000,    0.3176,    0.6196, 1]; % purple
color3 = [0.4660, 0.6740, 0.1880, 1]; % Green
color4 = [0.6235,    0.6000,    0.5961, 1]; % Red-Orange
color5 = [0.9290, 0.6940, 0.1250, 1];      % Blue

% Plot 1: Vertical sway angle
h1 = subplot(5,1,1);
xline(time(alpha), 'Color', [0, 0.447, 0.741, 0.7], 'LineStyle', '-.', 'LineWidth', 1.5); % Light gray dashed line
hold on;
xline(time(end-alpha+1), 'Color', [0, 0.447, 0.741, 0.7], 'LineStyle', '-.', 'LineWidth', 1.5); % Light gray dashed line
p4 = plot(time, wr(1:q:end), '--', 'LineWidth', 1.5, 'Color', color1);
p2 = plot(time, w_model(1:q:end), '-', 'LineWidth', 1.5, 'Color', color5);
p1 = plot(time, w_dd(1:q:end), '-', 'LineWidth', 1.5, 'Color', color4);
p3 = plot(time, w_opt(1:q:end), '-', 'LineWidth', 1.5, 'Color', color2);
scatter1 = scatter(time(1:15:end), w_opt(1:15*q:end), 'o', 'MarkerFaceColor', color2(1:3),'MarkerEdgeColor', color2(1:3), 'SizeData', 20);
scatter1.MarkerFaceAlpha = 0.5;
% Dummy line for the legend
dummyLine = plot(nan, nan, 'o-', 'LineWidth', 0.8, 'Color', color2);
ylabel('$\ddot{\theta}_4$ [rad/s$^2$]', 'FontSize', 14, 'Interpreter', 'latex');
legend([p4, dummyLine, p1, p2], {'Reference$\;$', 'Optimal trajectory$\;$ \newline', 'Nonparametric prediciton', 'Model-based prediciton'}, 'Location', 'northoutside', 'Orientation', 'horizontal', 'Interpreter', 'latex', 'FontSize', 14, 'Box', 'off');
hold off;
set(gca,'FontSize', 14)
set(gca, 'XTick', []); % Hide ticks
set(h1, 'Position', positions(1, :));
ylim([-0.1, 0.1])
box on;

% Plot 2: Horizontal sway angle
h2 = subplot(5,1,2);
p4 = plot(time, wr(2:q:end), '--', 'LineWidth', 1.5, 'Color', color1);
hold on;
p2 = plot(time, w_model(2:q:end), '-', 'LineWidth', 1.5, 'Color', color5);
p1 = plot(time, w_dd(2:q:end), '-', 'LineWidth', 1.5, 'Color', color4);
p3 = plot(time, w_opt(2:q:end), '-', 'LineWidth', 1.5, 'Color', color2);
scatter1 = scatter(time(1:15:end), w_opt(2:15*q:end), 'o', 'MarkerFaceColor', color2(1:3),'MarkerEdgeColor', color2(1:3), 'SizeData', 20);
scatter1.MarkerFaceAlpha = 0.5;

xline(time(alpha), 'Color', [0, 0.447, 0.741, 0.7], 'LineStyle', '-.', 'LineWidth', 1.5); % Light gray dashed line
xline(time(end-alpha+1), 'Color', [0, 0.447, 0.741, 0.7], 'LineStyle', '-.', 'LineWidth', 1.5); % Light gray dashed line
hold off;
ylabel('$\theta_1$ [rad]', 'FontSize', 14, 'Interpreter', 'latex');
set(gca, 'XTick', []); % Hide ticks
set(gca,'FontSize',14)
set(h2, 'Position', positions(2, :));
ylim([-0.01, 0.01])

% Plot 3: Another angle or parameter
h3 = subplot(5,1,3);
p4 = plot(time, wr(3:q:end), '--', 'LineWidth', 1.5, 'Color', color1);
hold on;
p2 = plot(time, w_model(3:q:end), '-', 'LineWidth', 1.5, 'Color', color5);
p1 = plot(time, w_dd(3:q:end), '-', 'LineWidth', 1.5, 'Color', color4);
p3 = plot(time, w_opt(3:q:end), '-', 'LineWidth', 1.5, 'Color', color2);
scatter1 = scatter(time(1:15:end), w_opt(3:15*q:end), 'o', 'MarkerFaceColor', color2(1:3),'MarkerEdgeColor', color2(1:3), 'SizeData', 20);
scatter1.MarkerFaceAlpha = 0.5;
xline(time(alpha), 'Color', [0, 0.447, 0.741, 0.7], 'LineStyle', '-.', 'LineWidth', 1.5); % Light gray dashed line
xline(time(end-alpha+1), 'Color', [0, 0.447, 0.741, 0.7], 'LineStyle', '-.', 'LineWidth', 1.5); % Light gray dashed line
hold off;
ylabel('$\theta_2$ [rad]', 'FontSize', 14, 'Interpreter', 'latex');
set(gca, 'XTick', []); % Hide ticks
set(gca,'FontSize',14)
set(h3, 'Position', positions(3, :));
ylim([-0.01, 0.01])

% Plot 4: Horizontal boom angle
h4 = subplot(5,1,4);
p4 = plot(time, wr(4:q:end), '--', 'LineWidth', 1.5, 'Color', color1);
hold on;
p2 = plot(time, w_model(4:q:end), '-', 'LineWidth', 1.5, 'Color', color5);
p1 = plot(time, w_dd(4:q:end), '-', 'LineWidth', 1.5, 'Color', color4);
p3 = plot(time, w_opt(4:q:end), '-', 'LineWidth', 1.5, 'Color', color2);
scatter1 = scatter(time(1:15:end), w_opt(4:15*q:end), 'o', 'MarkerFaceColor', color2(1:3),'MarkerEdgeColor',color2(1:3), 'SizeData', 20);
scatter1.MarkerFaceAlpha = 0.5;
xline(time(alpha), 'Color', [0, 0.447, 0.741, 0.7], 'LineStyle', '-.', 'LineWidth', 1.5); % Light gray dashed line
xline(time(end-alpha+1), 'Color', [0, 0.447, 0.741, 0.7], 'LineStyle', '-.', 'LineWidth', 1.5); % Light gray dashed line
hold off;
ylabel('$\theta_4$ [rad]', 'FontSize', 14, 'Interpreter', 'latex');
set(gca, 'XTick', []); % Hide ticks
set(gca,'FontSize',14)
set(h4, 'Position', positions(4, :));
ylim([wr(4)-0.05, wr(end-1)+0.05])

% Plot 5: Another parameter
h5 = subplot(5,1,5);
p4 = plot(time, wr(5:q:end), '--', 'LineWidth', 1.5, 'Color', color1);
hold on;
p2 = plot(time, w_model(5:q:end), '-', 'LineWidth', 1.5, 'Color', color5);
p1 = plot(time, w_dd(5:q:end), '-', 'LineWidth', 1.5, 'Color', color4);
p3 = plot(time, w_opt(5:q:end), '-', 'LineWidth', 1.5, 'Color', color2);
scatter1 = scatter(time(1:15:end), w_opt(5:15*q:end), 'o', 'MarkerFaceColor', color2(1:3),'MarkerEdgeColor', color2(1:3), 'SizeData', 20);
scatter1.MarkerFaceAlpha = 0.5
xline(time(alpha), 'Color', [0, 0.447, 0.741, 0.7], 'LineStyle', '-.', 'LineWidth', 1.5); % Light gray dashed line
xline(time(end-alpha+1), 'Color', [0, 0.447, 0.741, 0.7], 'LineStyle', '-.', 'LineWidth', 1.5); % Light gray dashed line
hold off;
xlabel('Time [s]', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$\dot{\theta}_4$ [rad/s]', 'FontSize', 14, 'Interpreter', 'latex');
set(gca,'FontSize',14)
set(h5, 'Position', positions(5, :));
ylim([-0.05, 0.15])

% Adjust figure dimensions
width = 22.5;  % centimeters
height = 27; % centimeters
set(fig,'Units','centimeters')
set(fig, 'Position', [0, 0, width, height]); % [left, bottom, width, height]
set(fig, 'PaperSize', [width, height]);

%%
function cost = run_traj_generation(params, H_w,  w_vec_given, dt, m, n, L, q, coefs)
    % Unpack mu and sigma from the structure
    mu = params.mu;
    sigma = params.sigma;
    lambda = params.lambda;
    
    % Call the original function
    [~, cost] = generate_trajectory(H_w, w_vec_given, dt, m, n, L, q, coefs, lambda, mu, sigma, 0, 0);
end

