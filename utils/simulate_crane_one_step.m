% Integrate the crane dynamics by one time step using ode45.
%
% The full state is x = [theta1; dtheta1; theta2; dtheta2; theta4; dtheta4]
%   theta1  - load sway angle in the radial direction     [rad]
%   dtheta1 - rate of theta1                              [rad/s]
%   theta2  - load sway angle in the tangential direction [rad]
%   dtheta2 - rate of theta2                              [rad/s]
%   theta4  - horizontal boom angle                       [rad]
%   dtheta4 - horizontal boom angular velocity            [rad/s]
%
% Inputs:
%   old_state - 6x1 state vector at the current step
%   newInput  - boom angular acceleration ddtheta4 [rad/s^2] (control input)
%   dt        - integration time step [s]
%   noise     - logical flag; if true, adds measurement noise to the output
%
% Outputs:
%   old_state - 6x1 updated state vector after one step
%   newY      - 4x1 measured output [theta1; theta2; theta4; dtheta4]
%
% Physical parameters (BehavioralCrane laboratory setup):
%   g   = 9.81   m/s^2    gravitational acceleration
%   th3 = 0.8108 rad      fixed vertical boom elevation angle
%   L   = 2.2511 m        effective boom arm length
%   l   = 1.0    m        rope length

function [old_state, newY] = simulate_crane_one_step(old_state, newInput, dt, noise)
    dynamics = @(~, x) crane_dynamics(x, newInput);
    [~, Xtemp] = ode45(dynamics, [0, dt], old_state, odeset('RelTol', 1e-8, 'AbsTol', 1e-10));
    old_state = Xtemp(end, :)';
    newY = [old_state(1); old_state(3); old_state(5); old_state(6)];
    if noise
        W    = diag([0.005, 0.005, 0.01, 0.01]);
        newY = newY + deg2rad(W * randn(4, 1));
    end
end

% Full 6-state crane dynamics (right-hand side for ode45).
%   x       - state vector [theta1; dtheta1; theta2; dtheta2; theta4; dtheta4]
%   ddtheta4 - control input (boom angular acceleration) [rad/s^2]
function x_dot = crane_dynamics(x, ddtheta4)
    g   = 9.81;    % gravitational acceleration [m/s^2]
    th3 = 0.8108;  % fixed vertical boom elevation angle [rad]
    L   = 2.2511;  % effective boom arm length [m]
    l   = 1.0;     % rope length [m]

    theta1  = x(1);
    dtheta1 = x(2);
    theta2  = x(3);
    dtheta2 = x(4);
    dtheta4 = x(6);  % x(5) = theta4, not needed explicitly below

    Vel = dtheta4;   % boom angular velocity  (dtheta4)
    Acc = ddtheta4;  % boom angular acceleration (control input)

    x_dot    = zeros(6, 1);
    x_dot(1) = dtheta1;
    x_dot(2) = ( -l^2*(-2*dtheta2*Vel - theta2*Acc - theta1*Vel^2 + theta1*dtheta1^2 + theta1*dtheta2^2) ...
                 -L*l*(-Vel^2*sin(th3)) ...
                 - g*l*theta1 ) / (l^2 + l^2*theta1^2);
    x_dot(3) = dtheta2;
    x_dot(4) = ( -L*l*(Acc*sin(th3)) ...
                 -l^2*(2*dtheta1*Vel + theta1*Acc + theta2*dtheta2^2 + dtheta1^2*theta2 - theta2*Vel^2) ...
                 - g*l*theta2 ) / (l^2 + l^2*theta2^2);
    x_dot(5) = dtheta4;
    x_dot(6) = ddtheta4;
end
