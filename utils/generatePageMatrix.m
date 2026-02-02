function [P, ud, yd, wd] = generatePageMatrix(numberOfTraj, L, q, dt, noise)
    P = [];
    maxAmplitude = abs(3/(30/pi*243*(-1/100))); % 2V changed to rad/s;
    T = 1200;
    old_state = [0; 0; 0; 0; pi/4; 0];
    time = 0:0.05:((T-1)*0.05);
    halt = zeros(100,1);

    ud = [0];  % Input data matrix
    yd = [old_state(1), old_state(3), old_state(5), old_state(6)];  % Output data matrix
    wd = [ud, yd]; % Data matrix
    k = 0;
    while (size(P,2) < numberOfTraj)
        k = k+1;
        old_state(5) = pi/4;
        dtheta4 = generate_random_input(T, maxAmplitude);
        u = calculateAcceleration(dtheta4, time);
        %u = [u; halt];
        sub_w = [];
        for i=1:length(u)
            newInput = u(i)+0.0008*(rand-0.5);
            [old_state, newY] = simulate_crane_one_step(old_state, newInput, dt, noise);
            if i <= T  || k==1
                % Append new data points
                ud = [ud; newInput];
                yd = [yd; newY'];
    
                new_w = [newInput, newY'];
                sub_w = [sub_w; new_w];
                wd = [wd; new_w]; 
            end
        end
        H_w = blkhank(sub_w, L);
        P = [P, H_w];
    end
    disp(['number of separate trajectories: ', num2str(k)]);
    end


    
   
%     figure;
%     subplot(2,1,1)
%     plot(time, dtheta4);
%     subplot(2,1,2)
%     plot(time, u);

    
        


