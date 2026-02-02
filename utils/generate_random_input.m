% generate random input from summing up 20 sine function, smoothing the
% function and taper it at the end.

function input = generate_random_input(numPoints, maxAmplitude)

    % Generate random values between -0.1 and 0.1
    rand_freq = (10*randn(1, 20) - 5);
    t = 0:0.05:(numPoints-1)*0.05;
    randomValues = 0;
    for i = 1:20
        randomValues = randomValues + sin(rand_freq(i)*t);
    end
    
    % Apply a two-pass moving average filter to smooth the trajectory
    windowSize = 70; % 70 % Adjust this value to control the smoothness
    smoothTrajectory = filtfilt(ones(1, windowSize)/windowSize, 1, randomValues);
    
    % Normalize the trajectory to ensure the maximum absolute value is 0.1
    maxValue = max(abs(smoothTrajectory));
    smoothTrajectory = smoothTrajectory * maxAmplitude / maxValue;

    % Apply tapering to the start and end of the trajectory
    taperLength = floor(numPoints * 0.15); % 10% of the total number of points
    taper = hann(taperLength*2);
    taperStart = taper(1:taperLength);    % First half of Hann window
    taperEnd = taper(taperLength+1:end);  % Second half of Hann window

    % Apply taper to the trajectory
    smoothTrajectory(1:taperLength) = smoothTrajectory(1:taperLength) .* taperStart';
    smoothTrajectory(end-taperLength+1:end) = smoothTrajectory(end-taperLength+1:end) .* taperEnd';

    input = smoothTrajectory';

%     % Plot the trajectory
%     figure;
%     plot(input);
%     title('Random Smooth Trajectory');
%     xlabel('Time');
%     ylabel('Amplitude');

end