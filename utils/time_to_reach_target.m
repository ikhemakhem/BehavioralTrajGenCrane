% Return the first time the boom stays close to the target long enough.
function t = time_to_reach_target(w, target, tol, duration, dt, q)
    % Extract the relevant entries from w
    theta4_vector = w(4:q:end);
    
    % Calculate the absolute difference from the target
    theta4_diff = abs(theta4_vector - target);
    
    % Find indices where the difference is within the tolerance
    within_tolerance = theta4_diff <= tol;
    
    % Count consecutive time steps within tolerance
    if any(within_tolerance)
        % Find continuous stretches where the condition is true
        %starts = find(diff(within_tolerance(1:end-1)) == 0 & within_tolerance) == 1);
        %ends = find(diff([within_tolerance, 0] == 0 & [0, within_tolerance]) == 1);
        transitions = diff([0; within_tolerance; 0] == 1);
        starts = find(transitions == 1);
        ends = find(transitions == -1) - 1;

        % Calculate lengths of 1-sequences
        lengths = ends - starts + 1;

        % Find indices of sequences that exceed the duration
        long_seq_indices = find(lengths*dt > duration);

         % Determine start position based on condition
        if ~isempty(long_seq_indices)
            % Return the start of the first long sequence
            t = (starts(long_seq_indices(1))-1)*dt;
        else
            % Return the start of the last sequence if no long sequence is found
            if ~isempty(starts)
                t = (starts(end)-1)*dt;
            else
                % Return NaN if there are no sequences of 1s
                t = inf;
            end
        end
    else
        t = inf;
    end
end
