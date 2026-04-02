% Score a generated trajectory using time, sway, smoothness, and overshoot.

function cost = quantify_function(w, target, T_target, dt, L, q)
    u = w(1:q:end);
    theta1 = w(2:q:end);
    theta2 = w(3:q:end);
    theta4 = w(4:q:end);
    %dtheta4dt = w(5:q:end);
    u_zero_padded = [zeros(2,1); u; zeros(2,1)];
    theta4_zero_padded = [theta4(1)*ones(2,1); theta4; theta4(end)*ones(2,1)];
    % compute derivative of input
    D = zeros(length(u), length(u_zero_padded));
    for i = 1:length(u)
        D(i,i) = 1/12;
        D(i,i+1) = -2/3;
        D(i,i+2) = 0;
        D(i,i+3) = 2/3;
        D(i,i+4) = -1/12;
    end
    
    dudt = (1/dt)*D*u_zero_padded;
    average_dudt = mean(abs(dudt));
    max_dudt = max(abs(dudt));

    dtheta4dt = (1/dt)*D*theta4_zero_padded;
    sum_dtheta4dt = sum(abs(dtheta4dt));
    max_dtheta4dt = max(abs(dtheta4dt));

    % measure oscilations from reference value
    max_oscillation = max(max(abs(theta1)),max(abs(theta2)));
    average_oscillation = (mean(abs(theta1)) + mean(abs(theta2)))/2;

    % Overshooting
    overshooting = sum(max(0, theta4 - target));

    % Normalization
    % Expected max or typical ranges
    % Normalization factors
    max_dudt_expected = 27.02;%10.2668;
    average_dudt_expected =  7.21;%
    sum_dtheta4dt_expected =  287.9;
    max_oscillation_expected = 0.070;
    average_oscillation_expected = 0.022;
    overshooting_expected = 3.51;

    % Normalize the terms
    norm_T_target = T_target / (L*dt);
    norm_max_dudt = max_dudt / max_dudt_expected;
    norm_average_dudt = average_dudt / average_dudt_expected;
    
    norm_sum_dtheta4dt = sum_dtheta4dt / sum_dtheta4dt_expected;
    norm_max_oscillation = max_oscillation / max_oscillation_expected;
    norm_average_oscillation = average_oscillation / average_oscillation_expected;
    norm_overshooting = overshooting / overshooting_expected;

    % Weights
    k1 = 3 ;
    k2 = 1;
    k3 = 3;
    k4 = 1/3;
    k5 = 1;

    % Calculate cost
    cost = k1 * norm_T_target + ...
           k2 * (2*norm_max_dudt + norm_average_dudt) + ...
           k3 * norm_sum_dtheta4dt + ...
           k4 * (2*norm_max_oscillation + norm_average_oscillation) + ...
           k5 * norm_overshooting;

    %%

    if isnan(cost) || (max_dtheta4dt > 5) 
        cost = 1000;
    end

    disp(['values: T_target, max_dudt, average_dudt, max_dtheta4dt, norm_sum_dtheta4dt, max_oscillation, average_oscillation, overshooting:  ', num2str(T_target),' , ', num2str(max_dudt),' , ', num2str(average_dudt),' , ', num2str(max_dtheta4dt),' , ', num2str(norm_sum_dtheta4dt),' , ', num2str(max_oscillation),' , ', num2str(average_oscillation),' , ', num2str(overshooting)]);
end
