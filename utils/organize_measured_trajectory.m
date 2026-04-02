% Reconstruct and filter one measured trajectory from experimental results.
function w_meas = organize_measured_trajectory(results_data)
    cutoffFreq = 8/20;
    [lpb, lpa] = butter(5, cutoffFreq, 'low');

    theta_3 = pi/2 - 0.76;
    rope_length = 1;
    theta_4_offset = 0;
    z_boom = 1.55;

    raw_vision_data = results_data.raw_vision_data;
    vision_data = results_data.vision_data;

    dtheta_4 = cellfun(@(c) c{2}, vision_data);
    ddtheta_4 = cellfun(@(c) c{3}, vision_data);
    t4 = cellfun(@(c) c{4}, raw_vision_data);
    u = cellfun(@(c) c{5}, raw_vision_data);
    v = cellfun(@(c) c{6}, raw_vision_data);

    numSamples = numel(u);
    theta_1_ = zeros(numSamples, 1);
    theta_2_ = zeros(numSamples, 1);
    theta_4_ = zeros(numSamples, 1);

    for idx = 1:numSamples
        [theta1_i, theta2_i, theta4_i] = pixels_to_world( ...
            u(idx), v(idx), t4(idx), theta_3, rope_length, theta_4_offset, z_boom);
        theta_1_(idx) = 1.7 * theta1_i;
        theta_2_(idx) = 1.7 * theta2_i;
        theta_4_(idx) = theta4_i;
    end

    theta_1_ = theta_1_ - mean(theta_1_);
    theta_2_ = theta_2_ - mean(theta_2_);

    theta_1 = filtfilt(lpb, lpa, theta_1_);
    theta_2 = filtfilt(lpb, lpa, theta_2_);
    theta_4 = linear_kalman(ddtheta_4, [theta_4_, dtheta_4]);
    theta_4 = theta_4(:,1);

    w_meas = [ddtheta_4, theta_1, theta_2, theta_4, dtheta_4]';
    w_meas = w_meas(:);
end
