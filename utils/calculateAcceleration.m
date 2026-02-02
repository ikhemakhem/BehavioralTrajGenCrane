% compute derivative of the velocity to compute the acceleration.
function acceleration = calculateAcceleration(velocity, time)
    % This function calculates the acceleration from a given velocity array using a fourth-order central difference method.
    % Input:
    %   velocity - A vector of velocities
    %   time     - A vector of times corresponding to the velocities
    % Output:
    %   acceleration - A vector of calculated accelerations
    
    n = length(velocity);            % Number of data points
    acceleration = zeros(n, 1);      % Preallocate the acceleration array

    % Check if time vector is of the same length as velocity vector
    if length(time) ~= n
        error('Time and velocity vectors must be of the same length');
    end

    % Calculate the acceleration using numerical differentiation
    for i = 1:n
        if i == 1 || i == 2
            % Second-order forward difference for the first two points
            h1 = time(i+1) - time(i);
            h2 = time(i+2) - time(i);
            acceleration(i) = (-velocity(i+2) + 4*velocity(i+1) - 3*velocity(i)) / (2*h1);
        elseif i == n || i == n-1
            % Second-order backward difference for the last two points
            h1 = time(i) - time(i-1);
            h2 = time(i) - time(i-2);
            acceleration(i) = (3*velocity(i) - 4*velocity(i-1) + velocity(i-2)) / (2*h1);
        else
            % Fourth-order central difference for the middle points
            h = (time(i+1) - time(i-1)) / 2;  % Assuming uniform spacing
            acceleration(i) = (-velocity(i+2) + 8*velocity(i+1) - 8*velocity(i-1) + velocity(i-2)) / (12*h);
        end
    end
end
