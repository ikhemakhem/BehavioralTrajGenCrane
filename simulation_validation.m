% simulation_validation.m
% Validate the nonparametric behavioral simulation framework on simulated data.
%
% Generates crane trajectories from the model, checks the Willems/Berberich
% rank conditions, and compares the nonparametric prediction against the
% model output for a random test input sequence.
%
% Author: Iskandar Khemakhem
% Year: 2024

%% ----------------------------- Setup -----------------------------------
scriptDir = fileparts(mfilename('fullpath'));
addpath(genpath(scriptDir));

%% ----------------------------- System parameters -----------------------
dt = 0.05;  % sampling time [s]
m  = 1;     % number of inputs
n  = 12;    % upper estimate of system order
L  = 300;   % trajectory length (time steps)
q  = 5;     % stacked sample size: [ddtheta4, theta1, theta2, theta4, dtheta4]

%% ----------------------------- Generate simulated data -----------------
numberOfTraj = 12 * q * L;
noise        = true;
[H_w, ud, yd, ~] = generatePageMatrix(numberOfTraj, L, q, dt, noise);

figure
subplot(2,3,1); plot(rad2deg(ud));
title('Simulated: input $\ddot{\theta}_4$', 'Interpreter','latex')

subplot(2,3,2); plot(rad2deg(yd(:,1)));
title('Simulated: output $\theta_1$', 'Interpreter','latex')

subplot(2,3,3); plot(rad2deg(yd(:,2)));
title('Simulated: output $\theta_2$', 'Interpreter','latex')

subplot(2,3,5); plot(rad2deg(yd(:,3)));
title('Simulated: output $\theta_4$', 'Interpreter','latex')

subplot(2,3,6); plot(rad2deg(yd(:,4)));
title('Simulated: output $\dot{\theta}_4$', 'Interpreter','latex')

% Verify rank conditions
disp(['More data than L?: ', num2str(length(ud) > L)]);
H_u_willems   = blkhank(ud, L + n);
H_u_berberich = blkhank(ud, L + n + 1);
disp(['Willems condition satisfied? ',   num2str(rank(H_u_willems)  >= m*L + n)]);
disp(['Berberich condition satisfied? ', num2str(rank(H_u_berberich) >= m*(L+n+1))]);

%% ----------------------------- Define test trajectory -----------------
init_cond     = 3;
init_cond_idx = 1:q*init_cond;
input_idx     = q*init_cond + 1 : q : q*L;
I_given       = [init_cond_idx, input_idx];

w_init     = repmat([0; 0; 0; pi/4; 0], init_cond, 1);
input_seq  = generate_random_input(L - init_cond, 0.15);

H_w_given    = H_w(I_given, :);
w_vec_given  = [w_init; input_seq];

disp(['Is given data consistent?: ', ...
      num2str(rank([w_vec_given, H_w_given]) == rank(H_w_given))]);

%% ----------------------------- Compare simulations --------------------
threshold = 0.01;
lambda    = 1.838974522562955e-04;

w_dd    = simulate_dd_system(H_w, I_given, w_vec_given, m, n, L, q, lambda, init_cond);
w_model = simulate_model(w_vec_given, dt, init_cond, q);

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
