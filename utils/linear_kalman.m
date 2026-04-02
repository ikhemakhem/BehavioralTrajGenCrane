% Estimate boom angle and velocity with a simple linear Kalman filter.
function x_est_store = linear_kalman(u_series, y_series)

    % Define the sampling time
    T = 0.05;  % Define your sampling time here, e.g., 0.1 seconds
    
    % Define the continuous system dynamics
    A = [0 1; 0 0];  % State transition matrix
    B = [0; 1];      % Control input matrix
    C = [1 0; 0 1];  % Observation matrix (observe both position and velocity)
    D = 0;           % Direct transmission matrix
    
    % Convert continuous system to discrete system
    sys_continuous = ss(A, B, C, D);
    sys_discrete = c2d(sys_continuous, T);
    [Ad, Bd, Cd, Dd] = ssdata(sys_discrete);
    
    % Define the covariance matrices
    Q = [0 0; 0.0 0];  % Process noise covariance
    R = [200 0; 0 0.001]; % Measurement noise covariance
    
    % Initialize the state estimate and covariance
    x_est = y_series(1,:)';  % Initial state estimate
    P = [0.1 0; 0 0.001];         % Initial estimate covariance
    
    % Preallocate arrays for storing results
    num_steps = length(u_series);
    x_est_store = zeros(num_steps, 2);  % Estimated state
    
    % Apply the Kalman filter to the real data
    for k = 1:num_steps
        % Control input for current step
        u = u_series(k);
        
        % Noisy measurement for current step
        y = y_series(k,:)';
        
        % Prediction step
        x_pred = Ad * x_est + Bd * u;
        P_pred = Ad * P * Ad' + Q;
        
        % Update step
        K = P_pred * Cd' / (Cd * P_pred * Cd' + R);  % Kalman gain
        x_est = x_pred + K * (y - Cd * x_pred);
        P = (eye(2) - K * Cd) * P_pred;
        
        % Store the estimated state
        x_est_store(k, :) = x_est';
    end
    
%     % Plot the results
%     figure;
%     subplot(2, 1, 1);
%     plot(1:num_steps, y_series(:,1), 'r-', 'LineWidth', 1);
%     hold on;
%     plot(1:num_steps, x_est_store(:, 1), 'b', 'LineWidth', 1);
%     legend('Noisy Observations theta_4', 'Estimated State theta_4');
%     title('State theta_4 Estimation');
%     xlabel('Time step');
%     ylabel('\theta_4');
% 
%     subplot(2, 1, 2);
%     plot(1:num_steps, y_series(:,2), 'r--', 'LineWidth', 1);
%     hold on;
%     plot(1:num_steps, x_est_store(:, 2), 'b', 'LineWidth', 1);
%     legend('Noisy Observations dtheta4', 'Estimated State dthtea4');
%     title('State dtheta4 Estimation');
%     xlabel('Time step');
%     ylabel('dtheta4');
    
end

