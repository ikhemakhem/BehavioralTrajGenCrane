% simulate crane dynamics and output the results in the same format of a
% trajectory in the behavioral settings
function w_dd = simulate_dd_system(H_w, I_given, w_vec_given, m, n, L, q, lambda, init_cond)
    H_w_given = H_w(I_given,:);  
        %% Data-driven interpolation
        % Weight matrix
        W = diag([1*ones(q,1);ones(length(I_given)-q,1)]);
        

%         lambda = 0.01

%         B = [eye(q*init_cond), zeros(q*init_cond,length(w_vec_given)-q*init_cond);zeros(5,length(w_vec_given)-5), eye(5)];
%         init = [w_vec_given(1:q*init_cond); w_vec_given(end-4:end)];

        seq_len = length(w_vec_given);
        B = [eye(q*init_cond), zeros(q*init_cond, seq_len-q*init_cond)];
        init = [w_vec_given(1:q*init_cond)];
        C = [zeros(seq_len-q*init_cond, q*init_cond), eye(seq_len-q*init_cond)];


%         lambda = 0;
        
        cvx_begin
        
            variable g(size(H_w_given, 2)); % Replace with the actual size of g
            minimize(1/2*(w_vec_given-H_w_given*g)'*W*(w_vec_given-H_w_given*g) + lambda * norm(g,1));
            subject to
                %ones(1, size(H_w_given, 2)) * g == 1
                abs(B*H_w_given*g - init) <= 0.00001*ones(length(init),1)
                abs(C*H_w_given*g - C*w_vec_given) <= 0.0001*ones(seq_len-q*init_cond,1)
        cvx_end

        %g = ((sqrtm(W)*H_w_given)\sqrtm(W))*w_vec_given;
        
        g_opt = g; % Optimal g
        lasso_val = cvx_optval; % Optimal value of the LASSO cost function
        
        % check sparsity of g_opt
%         if threshold > 0
%             disp(['rank of H_w after truncation: ', num2str(p)]);
%             disp(['r ', num2str(r)]);
%         end
        disp(['Will solution be unique?: ' num2str(rank(H_w_given) == rank(H_w))]);    % Test uniqueness
        disp(['sum of g ', num2str(ones(1, size(H_w_given, 2)) * g_opt)]);
        disp(['Nuclear norm upper value ', num2str(L*m+n)]);
        disp(['Sparsity of optimal solution ', num2str(sum(abs(g_opt) > (1e-4)*ones(size(g_opt))))]);

        fprintf('\n');
        w_dd = H_w*g_opt;
end