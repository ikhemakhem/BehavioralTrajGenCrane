%% Validate framework on Simulation
dt = 0.05; % time step
m = 1; 
n = 12; % upper esimation
L = 300;
q = 5; % input + output = size(w(t))

%% 1) Generate data
numberOfTraj = 12*q*L;
noise = true;
[H_w, ud, yd, wd] = generatePageMatrix(numberOfTraj, L, q, dt, noise);

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

% Verify
disp(['Do we have more data than L?: ', num2str(length(ud)> L)]);
% Test input condition
H_u_willems = blkhank(ud, L+n);
disp(['Is Willems Lemma condition satisfied? ', num2str(rank(H_u_willems) >= m*L+n)]);
H_u_berberich = blkhank(ud, L+n+1);
disp(['Is Berberich condition satisfied? ', num2str(rank(H_u_berberich) >= m*(L+n+1))]);

%% 2) Define given data
% init condition and input sequence
% make initial conditions long enough T_init > lag to to specify systems
% initital conditions
init_cond = 3;
init_cond_idx = 1:q*init_cond;
input_idx = q*init_cond+1:q:q*L;

I_given = [init_cond_idx, input_idx];
w_init = repmat([0; 0; 0; pi/4; 0], init_cond, 1);

% t = linspace(0,L*dt,L-init_cond);
% input_seq = 0.03*(sin(2 * pi*t')+ sin(pi*t')) + 0.001 * (randn(L-init_cond,1)-0.5); % Sine wave with noise
% input_seq = input_seq - input_seq(1)*ones(L-init_cond,1);

input_seq = generate_random_input(L-init_cond, 0.15);

H_w_given = H_w(I_given,:);

w_vec_given = [w_init; input_seq];

% Test consistency
disp(['Is given data consistent?: ', num2str(rank([w_vec_given, H_w_given]) == rank(H_w_given))]);

%% Compare simulations
threshold = 0.01;
lambda = 1.838974522562955e-04;
%lambda = 3e-04;
w_dd = simulate_dd_system(H_w, I_given, w_vec_given, m, n, L, q, threshold, lambda);
w_model = simulate_model(w_vec_given, dt, init_cond, q);

figure;
subplot(2,3,1);
plot(rad2deg(w_dd(1:q:end)));
hold on;
plot(rad2deg(w_model(1:q:end)));
title('Solution input ddtheta4');
legend('data driven', 'model');
subplot(2,3,2);
plot(rad2deg(w_dd(2:q:end)));
hold on;
plot(rad2deg(w_model(2:q:end)));
title('Solution output theta1');
legend('data driven', 'model');
subplot(2,3,3);
plot(rad2deg(w_dd(3:q:end)));
hold on;
plot(rad2deg(w_model(3:q:end)));
title('Solution output theta2');
legend('data driven', 'model');
subplot(2,3,5);
plot(rad2deg(w_dd(4:q:end)));
hold on;
plot(rad2deg(w_model(4:q:end)));
title('Solution output theta4');
legend('data driven', 'model');
subplot(2,3,6);
plot(rad2deg(w_dd(5:q:end)));
hold on;
plot(rad2deg(w_model(5:q:end)));
title('Solution output dtheta4');
legend('data driven', 'model');
hold off;



