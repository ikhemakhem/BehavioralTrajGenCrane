% Recover a full trajectory from partial observations using the behavioral model.
%
% Solves a weighted LASSO problem:
%   min  0.5*(w_given - H_given*g)'*W*(w_given - H_given*g) + lambda*||g||_1
%   s.t. initial-condition entries of H_w*g match w_vec_given to tolerance
%        remaining entries match the prescribed inputs to tolerance
%
% Requires CVX (http://cvxr.com/cvx).
%
% Inputs:
%   H_w         - (q*L x N) data matrix (block Hankel / Page matrix)
%   I_given     - indices of the known entries in the full trajectory vector
%   w_vec_given - values at the known indices (initial conditions + input sequence)
%   m           - number of inputs
%   n           - estimated system order
%   L           - trajectory length (time steps)
%   q           - stacked sample size (inputs + outputs)
%   lambda      - regularisation weight on ||g||_1
%   init_cond   - number of initial-condition time steps
%
% Output:
%   w_dd - (q*L x 1) reconstructed full trajectory

function w_dd = simulate_dd_system(H_w, I_given, w_vec_given, m, n, L, q, lambda, init_cond)
    H_w_given = H_w(I_given, :);

    seq_len = length(w_vec_given);
    W = diag([ones(q, 1); ones(length(I_given) - q, 1)]);
    B = [eye(q*init_cond),  zeros(q*init_cond, seq_len - q*init_cond)];
    C = [zeros(seq_len - q*init_cond, q*init_cond),  eye(seq_len - q*init_cond)];
    init = w_vec_given(1:q*init_cond);

    cvx_begin quiet
        variable g(size(H_w_given, 2));
        minimize( 1/2*(w_vec_given - H_w_given*g)' * W * (w_vec_given - H_w_given*g) ...
                  + lambda * norm(g, 1) );
        subject to
            abs(B * H_w_given * g - init)          <= 1e-5 * ones(length(init), 1)           % CVX constraint
            abs(C * H_w_given * g - C*w_vec_given) <= 1e-4 * ones(seq_len - q*init_cond, 1)  % CVX constraint
    cvx_end

    disp(['Will solution be unique?: ', num2str(rank(H_w_given) == rank(H_w))]);
    disp(['Sum of g: ',                 num2str(sum(g))]);
    disp(['Nuclear norm upper bound: ', num2str(L*m + n)]);
    disp(['Sparsity of g (|g|>1e-4): ', num2str(sum(abs(g) > 1e-4))]);
    fprintf('\n');

    w_dd = H_w * g;
end
