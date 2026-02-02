% Generate an optimal trajectory by solving a CVX optimization problem.
% Optionally plot the resulting signals and print diagnostic information.
%
% The decision variable is the coefficient vector g. The trajectory is then
% reconstructed as w_opt = H_w * g. The objective combines:
%   (i)   weighted data-consistency with the given entries (least-squares),
%   (ii)  an l1 penalty on g to promote sparsity (LASSO),
%   (iii) a smoothness/finite-difference penalty via H_tt (time-shifted difference),
%   (iv)  a weighted tracking term towards a reference wr.
% Constraints enforce boundary conditions and inequality limits.
%
% Inputs:
%   H_w                 - Data/Hankel-like matrix to reconstruct the full trajectory.
%   w_vec_given         - Vector of known/observed entries (with missing entries implicit).
%   dt                  - Sampling time [s].
%   m, n, L, q           - Problem dimensions (used for diagnostics and helper functions).
%   coefs               - Struct of precomputed matrices/vectors:
%                         .H_tt      : block-difference operator (with boundary padding)
%                         .H_w_given : submatrix mapping g to the known entries
%                         .wr        : reference trajectory (same space as H_w*g)
%                         .W         : weighting matrix for given-entry fit
%                         .B, .c     : equality constraint matrices/vectors (start/target)
%                         .M, .d     : inequality constraint matrices/vectors (limits)
%                         .Q         : weighting matrix for reference tracking term
%   lambda              - Weight on ||g||_1 (sparsity).
%   mu                  - Weight on smoothness term ||H_tt*[1;1;1;g;1;1;1]||^2.
%   sigma               - Weight on reference tracking term (normalized by ||H_w||_F).
%   show_fig            - If true, plot key signals of w_opt.
%   show_final_result   - If true, print rank/sparsity/consistency diagnostics.
%
% Outputs:
%   w_opt               - Optimal reconstructed trajectory (w_opt = H_w*g_opt).
%   cost                - Scalar cost used by the caller (currently set to cost_1).
function [w_opt, cost] = generate_trajectory(H_w,  w_vec_given, dt, m, n, L, q, coefs, lambda, mu, sigma, show_fig, show_final_result)

    % Display current regularization weights (useful for sweeps).
    disp(['lambda = ', num2str(lambda), ', mu = ', num2str(mu), ' sigma = ', num2str(sigma)]);
   
    % Unpack precomputed coefficient matrices/vectors from struct for readability.
    H_tt      = coefs.H_tt;      % block-difference operator (handles time-shift + boundaries)
    H_w_given = coefs.H_w_given; % mapping from g to the known entries in w_vec_given
    wr        = coefs.wr;        % reference trajectory
    W         = coefs.W;         % weights for fitting the given entries
    B         = coefs.B;         % equality constraint matrix (boundary conditions)
    c         = coefs.c;         % equality constraint vector
    M         = coefs.M;         % inequality constraint matrix (limits)
    d         = coefs.d;         % inequality constraint bounds
    Q         = coefs.Q;         % weights for reference tracking term
    
    % Solve convex optimization problem with CVX.
    % Decision variable: g (coefficients selecting columns of H_w / H_w_given).
    cvx_begin
        variable g(size(H_w_given, 2));

        % Objective:
        % 1) Weighted least squares fit to known entries:
        %    0.5 * || W^(1/2) * (w_vec_given - H_w_given*g) ||_2^2
        % 2) Sparsity promotion:
        %    lambda * ||g||_1
        % 3) Smoothness (time-difference) penalty on the stacked vector with padded boundaries:
        %    mu * || H_tt * [1;1;1; g; 1;1;1] ||_2^2
        % 4) Reference tracking term (scaled by Frobenius norm of H_w):
        %    sigma * (1/||H_w||_F) * || Q^(1/2) * (H_w*g - wr) ||_2^2
        minimize( ...
            1/2*(w_vec_given-H_w_given*g)'*W*(w_vec_given-H_w_given*g) + ...
            lambda*norm(g,1) + ...
            mu*(H_tt*[1;1;1;g;1;1;1])'*(H_tt*[1;1;1;g;1;1;1]) + ...
            sigma*(1/norm(H_w,'fro'))*(H_w*g-wr)'*Q*(H_w*g-wr) ...
        );

        % Constraints:
        % - Enforce start/target (boundary) conditions on the known-entry reconstruction.
        % - Enforce inequality limits (e.g., angle/velocity/acceleration bounds) on full trajectory.
        subject to
            B*H_w_given*g == c
            M*H_w*g <= d
    cvx_end
    
    % Extract optimal solution and objective value reported by CVX.
    g_opt     = g;
    lasso_val = cvx_optval; %#ok<NASGU> % stored for debugging if needed
    
    % Diagnostics: sparsity level of g_opt.
    disp(['Nuclear norm upper value ', num2str(L*m+n)]);
    disp(['Sparsity of optimal solution ', num2str(sum(abs(g_opt) > (1e-4)*ones(size(g_opt))))]);
    
    % Reconstruct optimal trajectory in original signal space.
    w_opt = H_w*g_opt;

    % Compute time to reach target based on the reference terminal value.
    target   = wr(end-1);
    T_target = time_to_reach_target(w_opt, target, deg2rad(2), 5, dt, q);

    % Compare reconstructed trajectory to a simulated model response.
    % Uses a short initialization window, then samples every q-th element afterwards.
    init_cond           = 20;
    I_given_simulation  = [[1:q*init_cond]'; [q*init_cond+1:q:length(w_opt)]'];
    w_given_simulation  = w_opt(I_given_simulation);
    w_model             = simulate_model(w_given_simulation, dt, init_cond, q);

    % Primary cost term (problem-specific metric).
    cost_1 = quantify_function(w_opt, target, T_target, dt, L, q);

    % Secondary cost term: discrepancy to simulation (with selected elements excluded).
    w_opt_filtered   = exclude_5q_plus_2_elements(w_opt);
    w_model_filtered = exclude_5q_plus_2_elements(w_model);
    cost_2           = sqrt((w_opt_filtered - w_model_filtered)' * (w_opt_filtered - w_model_filtered)); %#ok<NASGU>

    % Combined cost (currently overridden below).
    cost = cost_1 + 2*cost_2;
    cost = cost_1;   % NOTE: final output cost currently ignores cost_2.
    
    % Optional plotting of key components of the trajectory.
    if show_fig
        theta_limit = d(1);
        vel_limit   = d(5);

        figure;

        % Input ddtheta4 (assumed first channel).
        subplot(2,3,1);
        plot(rad2deg(w_opt(1:q:end)));
        xlabel('Time steps ($\Delta t = 0.05 \, s$)', 'Interpreter', 'latex');
        ylabel('Acceleration ($\mathrm{deg/s}^2$)', 'Interpreter', 'latex');
        title('Optimal solution input ddtheta4');

        % Output theta1 (assumed second channel) with bounds.
        subplot(2,3,2);
        plot(rad2deg(w_opt(2:q:end)));
        hold on;
        yline(rad2deg(theta_limit), 'r--');
        yline(rad2deg(-theta_limit), 'r--');
        ylim([-3,3]);
        hold off;
        xlabel('Time steps ($\Delta t = 0.05 \, s$)', 'Interpreter', 'latex');
        ylabel('Load sway angle in radial direction ($\mathrm{deg}$)', 'Interpreter', 'latex');
        title('Optimal solution output theta1');

        % Output theta2 (assumed third channel) with bounds.
        subplot(2,3,3);
        plot(rad2deg(w_opt(3:q:end)));
        hold on;
        yline(rad2deg(theta_limit), 'r--');
        yline(rad2deg(-theta_limit), 'r--');
        ylim([-3,3]);
        hold off;
        xlabel('Time steps ($\Delta t = 0.05 \, s$)', 'Interpreter', 'latex');
        ylabel('Load away angle in tangential direction ($\mathrm{deg}$)', 'Interpreter', 'latex');
        title('Optimal solution output theta2');

        % Output theta4 (assumed fourth channel) and mark time-to-target.
        subplot(2,3,5);
        plot(rad2deg(w_opt(4:q:end)));
        hold on;
        x_line = T_target/dt;
        y_min  = min(rad2deg(w_opt(4:q:end)));
        y_max  = max(rad2deg(w_opt(4:q:end)));
        plot([x_line x_line], [y_min y_max], 'r--');
        hold off;
        xlabel('Time steps ($\Delta t = 0.05 \, s$)', 'Interpreter', 'latex');
        ylabel('Horizontal boom velocity ($\mathrm{deg/s}$)', 'Interpreter', 'latex');
        title('final result output theta4');

        % Output dtheta4 (assumed fifth channel) with bounds.
        subplot(2,3,6);
        plot(rad2deg(w_opt(5:q:end)));
        hold on;
        yline(rad2deg(vel_limit), 'r--');
        yline(rad2deg(-vel_limit), 'r--');
        hold off;
        title('Particular solution output dtheta4');
    end

    % Print final scalar cost.
    disp(['cost = ', num2str(cost)]);

    % Optional detailed diagnostics about rank/consistency/uniqueness.
    if show_final_result
        disp(['Will solution be unique?: ' num2str(rank(H_w_given) == rank(H_w))]);
        disp(['sum of g ', num2str(ones(1, size(H_w_given, 2)) * g_opt)]);
        disp(['Nuclear norm upper value ', num2str(L*m+n)]);
        disp(['Sparsity of optimal solution ', num2str(sum(abs(g_opt) > (1e-4)*ones(size(g_opt))))]);
        disp(['Is optimal trajectory consistent?: ', num2str(rank([w_opt, H_w]) == rank(H_w))]);
    end
end

% Exclude every element whose index i satisfies mod(i,5)==2 by setting it to zero.
% (Despite the name, this excludes "2 mod 5" entries; q is not used here.)
function filtered_array = exclude_5q_plus_2_elements(array)
    indices_to_exclude = [];
    for i = 1:length(array)
        if mod(i, 5) == 2
            indices_to_exclude = [indices_to_exclude, i]; %#ok<AGROW>
        end
    end
    filtered_array = array;
    filtered_array(indices_to_exclude) = 0;
end