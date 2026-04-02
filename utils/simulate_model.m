% Simulate the crane model from initial conditions and a given input sequence.
function w_model = simulate_model(w_vec_given, dt, init_cond, q)

w_init = w_vec_given(1:q*init_cond);
input_seq = w_vec_given(q*init_cond+1:end);



%old_state = [w_init(2); 0; w_init(3); 0; w_init(4); w_init(5)];
old_state = [w_init(2); 0; w_init(3); 0; w_init(4); w_init(5)];
input = w_init(1:q:q*init_cond);
input = [input; input_seq];
w_model = [];
   for i=1:length(input)
        newInput = input(i);
        [old_state, newY] = simulate_crane_one_step(old_state, newInput, dt, 0);
        w_model = [w_model; newInput; newY(1); newY(2); newY(3); newY(4)];
   end
end
