% grid search for nonparametric simulation
%%
    dt = 0.05; % time ste
    m = 1; 
    n = 6;
    L = 500;
    q = 5;     % stacked sample size; includes ddtheta4 plus theta1, theta2, theta4, dtheta4
%%
    [H_w, ud, yd, wd] = pageMatrixData(L);

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

    %% load input sequences
    init_cond = 20;
    input_sequences = {};
    for i = 1:15
        input_pattern = sprintf('input_seq_sim_500_%d.mat', i);
        input_data = load(input_pattern);
        input_sequences{i} = input_data.acceleration(init_cond+1:end);
    end
    % Create a new figure
    figure;
    % Plot each sequence in a subplot
    for i = 1:15
        subplot(5, 3, i); % Arrange in 5 rows and 2 columns
        plot(input_sequences{i});
        title(sprintf('Acceleration Profile %d', i));
        xlabel('Time (s)');
        ylabel('Acceleration (m/s^2)');
        grid on;
    end
    
    %% load measured trajectories and reorganize them
    measured_trajectories = cell(15,1);

    for i=1:15
        results_pattern = sprintf('results_seq_sim_500_%d.mat', i);
        results_data = load(results_pattern);

        measured_trajectories{i} = organize_measured_trajectory(results_data);
    end
    %%
    % Define parameter ranges
    threshold_values = custom_spaced_values(5e-6, 0.005, 5, 4); %15  % 0.02 and 0.05 result to fail 0.006 and 0.01 result to bad results.
    lambda_values = custom_spaced_values(1e-8, 1e-3, 5, 10); %15
    p_values = 2500:1500:8000;            
    
    % Unique threshold and p values for precomputation
    threshold_unique = unique(threshold_values);
    p_unique = unique(p_values);
    
    % Precompute QR and then SVD components for each unique p and threshold
    svd_store = cell(numel(p_unique), numel(threshold_unique));
    for j = 1:numel(p_unique)
        p = p_unique(j);
        
        % Perform QR decomposition on H_w
        [~,~,P] = qr(H_w, 'vector');
        P = sort(P(1:p));  % Select the first p columns in sorted order
        H_w_trunc = H_w(:, P);
        
        for i = 1:numel(threshold_unique)
            threshold = threshold_unique(i);
            
            % Perform SVD on the truncated H_w
            [U, S, V] = svd(H_w_trunc);
            S_vec = diag(S);
            r = find(S_vec >= threshold, 1, 'last'); % Last significant singular value
            
            if isempty(r)
                r = size(U, 2); % Use all vectors if threshold is too small
            end

            H_w_trunc = U(:, 1:r)*S(1:r, 1:r)*V(:, 1:r)';

            % Store truncated U, S, V
            svd_store{j, i} = H_w_trunc;
        end
    end
%%    
    % Generate grid of all combinations
    [ThresholdGrid, LambdaGrid, PGrid] = ndgrid(threshold_values, lambda_values, p_values);
    
    % Flatten grids to enable vectorized evaluation
    threshold_vector = ThresholdGrid(:);
    lambda_vector = LambdaGrid(:);
    p_vector = PGrid(:);
    
    % Map indices for threshold and p in SVD store
    [~, threshold_idx] = ismember(threshold_vector, threshold_unique);
    [~, p_idx] = ismember(p_vector, p_unique);
%%
    % Define objective function evaluations
    objective_values = arrayfun(@(p_idx, t_idx, l) objective_function(svd_store{p_idx, t_idx}, l, ...
        input_sequences, measured_trajectories, m, n, L, q, init_cond), ...
        p_idx, threshold_idx, lambda_vector);
    
    % Find the index of the minimum objective value
    [minValue, idx] = min(objective_values);
    
    % Extract the best parameters
    best_params = struct('threshold', threshold_vector(idx), 'lambda', lambda_vector(idx), 'p', p_vector(idx));
    
    % Display the best parameters
    disp('Best Parameters:');
    disp(best_params);
    
    % Save the results
    save('sim_grid_search_real_data_500.mat', 'best_params', 'objective_values', 'threshold_values', 'lambda_values', 'p_values', 'threshold_vector', 'lambda_vector', 'p_vector', 'input_sequences', 'init_cond');

    %%
    % Plot the objective values
    figure;
    plot(objective_values, 'b-', 'DisplayName', 'Objective Values');
    hold on;
    
    % Identify and plot NaN values
    nan_indices = isnan(objective_values);
    plot(find(nan_indices), zeros(length(find(nan_indices)),1), 'rx', 'MarkerSize', 10, 'DisplayName', 'NaN Values');
    
    % Add labels and legend
    xlabel('Index');
    ylabel('Objective Value');
    title('Objective Values with NaN marked');
    legend show;
    hold off;
    %% Compare simulations with optimized parameters
    threshold = best_params.threshold;
    lambda = best_params.lambda;
    p = 2500; %best_params.p;
    
    % Perform QR decomposition on H_w
    [~,~,P] = qr(H_w, 'vector');
    P = sort(P(1:p));  % Select the first p columns in sorted order
    H_w_trunc = H_w(:, P);
    % Perform SVD on the truncated H_w
    [U, S, V] = svd(H_w_trunc);
    S_vec = diag(S);
    r = find(S_vec >= threshold, 1, 'last'); % Last significant singular value
    if isempty(r)
        r = size(U, 2); % Use all vectors if threshold is too small
    end
    H_w_trunc = U(:, 1:r)*S(1:r, 1:r)*V(:, 1:r)';
%%
    init_cond_idx = 1:q*init_cond;
    input_idx = q*init_cond+5:q:q*L;

    I_given = [init_cond_idx, input_idx];

    for i = 1:20

        input_pattern = sprintf('input_seq_sim_500_%d.mat', i);
        input_data = load(input_pattern);
        input_seq = input_data.velocity(init_cond+1:end);

        results_pattern = sprintf('results_seq_sim_500_%d.mat', i);
        results_data = load(results_pattern);
        w_meas = organize_measured_trajectory(results_data);

        w_vec_given = [w_meas(1:q*init_cond); input_seq];
        w_dd = simulate_dd_system(H_w_trunc, I_given, w_vec_given, m, n, L, q, lambda, init_cond);

        if threshold > 0
            disp(['rank of H_w after truncation: ', num2str(p)]);
            disp(['r ', num2str(r)]);
        end
    
        figure;
        subplot(2,3,1);
        plot(rad2deg(w_dd(1:q:end)));
        hold on;
        plot(rad2deg(w_meas(1:q:end)));
        xlabel('Time steps ($\Delta t = 0.05 \, s$)', 'Interpreter', 'latex');
        ylabel('Acceleration ($\mathrm{deg/s}^2$)', 'Interpreter', 'latex');
        title('Solution input ddtheta4');
        legend('data driven', 'measurements');
        subplot(2,3,2);
        plot(rad2deg(w_dd(2:q:end)));
        hold on;
        plot(rad2deg(w_meas(2:q:end)));
        xlabel('Time steps ($\Delta t = 0.05 \, s$)', 'Interpreter', 'latex');
        ylabel('Load sway angle in radial direction ($\mathrm{deg}$)', 'Interpreter', 'latex');
        title('Solution output theta1');
        legend('data driven', 'measurements');
        subplot(2,3,3);
        plot(rad2deg(w_dd(3:q:end)));
        hold on;
        plot(rad2deg(w_meas(3:q:end)));
        xlabel('Time steps ($\Delta t = 0.05 \, s$)', 'Interpreter', 'latex');
        ylabel('Load away angle in tangential direction ($\mathrm{deg}$)', 'Interpreter', 'latex');
        title('Solution output theta2');
        legend('data driven', 'measurements');
        subplot(2,3,5);
        plot(rad2deg(w_dd(4:q:end)));
        hold on;
        plot(rad2deg(w_meas(4:q:end)));
        xlabel('Time steps ($\Delta t = 0.05 \, s$)', 'Interpreter', 'latex');
        ylabel('Horizontal boom angle ($\mathrm{deg}$)', 'Interpreter', 'latex');
        title('Solution output theta4');
        legend('data driven', 'measurements');
        subplot(2,3,6);
        plot(rad2deg(w_dd(5:q:end)));
        hold on;
        plot(rad2deg(w_meas(5:q:end)));
        xlabel('Time steps ($\Delta t = 0.05 \, s$)', 'Interpreter', 'latex');
        ylabel('Horizontal boom velocity ($\mathrm{deg/s}$)', 'Interpreter', 'latex');
        title('Solution output dtheta4');
        legend('data driven', 'measurements');
        hold off;

    end
    %%     Define the objective function for Bayesian optimization
    % Evaluate one simulation hyperparameter tuple on measured trajectories.
    function cost = objective_function(svd_data, lambda, input_sequences, measured_trajectories, m, n, L, q, init_cond)
        %init_cond = params.init_cond; % Assuming init_cond is fixed as 3
%         threshold = params.threshold;
%         lambda = params.lambda;
%         p = params.p;

        total_cost = 0;

        init_cond_idx = 1:q*init_cond;
        input_idx = q*init_cond+1:q:q*L;
        
        I_given = [init_cond_idx, input_idx];
               
        % Run the simulation 15 times with different random inputs
        for i = 1:15
            input_seq = input_sequences{i};
            input_seq = input_seq(1:L-init_cond);
            measured_traj = measured_trajectories{i};
            
            w_vec_given = [measured_traj(1:q*init_cond); input_seq];

            w_dd = simulate_dd_system(svd_data, I_given, w_vec_given, m, n, L, q, lambda, init_cond);
            total_cost = total_cost + sqrt((w_dd - measured_traj)'*(w_dd - measured_traj));
        end

        cost = total_cost;
    end
