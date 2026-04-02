%%
%grid search for optimal trajectory generation 
%%
% Datd riven trajectory Generation
dt = 0.05; % time step
m = 1;
n = 6;
L = 500;
q = 5; % stacked sample size; includes ddtheta4 plus theta1, theta2, theta4, dtheta4

%% 1) Generate data
[H_w, ud, yd, wd] = pageMatrixFromData(L);

%%
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
y_start = [zeros(alpha,1), zeros(alpha,1), (3*pi/8)*ones(alpha,1), zeros(alpha,1)];
w_start = [u_start, y_start];
w_start = reshape(w_start', 1, q*alpha)';

u_end = zeros(alpha,1);
y_end = [zeros(alpha,1), zeros(alpha,1), (pi/2)*ones(alpha,1),zeros(alpha,1)];
w_end = [u_end, y_end];
w_end = reshape(w_end', 1, q*alpha)';

H_w_given = H_w(I_given,:);

w_vec_given = [w_start; w_end];

% Test consistency
disp(['Is given data consistent?: ', num2str(rank([w_vec_given, H_w_given]) == rank(H_w_given))]);

%% Load the results from simulation grid search
load('sim_grid_search_real_data_500.mat');
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
vel_limit = abs(4/(30/pi*243*(-1/100)));
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
mu_values = custom_spaced_values(1e-1, 50, 5, 6);   %1e-7, 100, 15, 6
sigma_values = custom_spaced_values(1e-1, 50, 5, 7); % 1e-2, 100, 7, 7
lambda_values = custom_spaced_values(1e-7, 1, 5, 7);
% mu_values = [1e-7, 1e-4, 5e-2,  1, 10, 25, 100];
% sigma_values = [1e-7,  5e-5, 5e-2,  1, 10, 25, 100];
% lambda_values = [1e-7, 5e-5, 5e-2, 1, 10, 25, 100];

% Generate grid of all combinations
[muGrid, sigmaGrid, LambdaGrid] = ndgrid(mu_values, sigma_values, lambda_values);

% Flatten grids to enable vectorized evaluation
mu_vector = muGrid(:);
sigma_vector = sigmaGrid(:);
lambda_vector = LambdaGrid(:);

%%
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
% load("opt_traj_1.mat");
% results_data = load('results_opt_traj_1.mat');
% 
% w_model = traj;
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
save('opt_traj_500_5.mat','w_opt', 'w_model', 'w_dd');
velocity = w_opt(5:q:end);
acceleration = w_opt(1:q:end);
position = w_opt(4:q:end);
save('opt_traj_500_5_input.mat','position', 'velocity', 'acceleration');

%%
% Evaluate one trajectory-generation hyperparameter tuple.
function cost = run_traj_generation(params, H_w,  w_vec_given, dt, m, n, L, q, coefs)
    % Unpack mu and sigma from the structure
    mu = params.mu;
    sigma = params.sigma;
    lambda = params.lambda;
    
    % Call the original function
    [~, cost] = generate_trajectory(H_w, w_vec_given, dt, m, n, L, q, coefs, lambda, mu, sigma, 0, 0);
end

