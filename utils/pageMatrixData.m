% Build a data matrix from raw experimental files in a local data folder.
function [P, ud, yd, wd] = pageMatrixData(L)
    cd('C:\Users\ilkhe\source\repos\Matlab\Data_Generation\data');
    valid_seq = [1,2,3,5,6,8,9,10,11,12,13,14,15,17,18,22,23,29,21,31,32,33,34,35,36,37,38,39,40];
    % 
    folder = 'test3';  % Specify the folder name
    P = [];
    ud = [];
    yd = [];
    wd = [];
    for i = valid_seq
        % Create the pattern for input_seq_i.mat
        input_pattern = sprintf('input_seq_%d.mat', i);
        input_file = fullfile(folder, input_pattern);
        
        % Load the input file if it exists
        if exist(input_file, 'file')
            input_data = load(input_file, 'velocity', 'acceleration');
            fprintf('Loaded: %s\n', input_file);
        else
            fprintf('Input file not found: %s\n', input_file);
            continue;
        end
        
        % Create the pattern for results_seq_i_x.mat with wildcard for x
        results_pattern = sprintf('results_seq_%d.mat', i);
        results_files = dir(fullfile(folder, results_pattern));
        
        % Check if any results files are found
        if ~isempty(results_files)
            for k = 1:length(results_files)
                results_file = fullfile(folder, results_files(k).name);
                results_data = load(results_file);
                fprintf('Loaded: %s\n', results_file);
            end
        else
            fprintf('Results files not found for sequence %d\n', i);
            continue;
        end
        
        velocity = input_data.velocity;
        T = length(velocity);
        time = 0:0.05:((T-1)*0.05);
        acceleration = calculateAcceleration(velocity, time);
        
        raw_vision_data = results_data.raw_vision_data;
        vision_data = results_data.vision_data;
        dtheta_4 = cellfun(@(c) c{2}, vision_data);
        ddtheta_4 = cellfun(@(c) c{3}, vision_data);
        t4 = cellfun(@(c) c{4}, raw_vision_data);
        u = cellfun(@(c) c{5}, raw_vision_data);
        v = cellfun(@(c) c{6}, raw_vision_data);
        
        theta_3 = pi/2 - 0.76;
        rope_length = 1; 
        theta_4_offset = 0;
        z_boom = 1.55;
        theta_4_ = [];
        theta_1_ = [];
        theta_2_ = [];
        for i = 1:length(u)
            [theta1_i, theta2_i, theta4_i] = pixels_to_world(u(i), v(i), t4(i), theta_3, rope_length, theta_4_offset, z_boom);
            theta_1_ = [theta_1_; 1.7*theta1_i];
            theta_2_ = [theta_2_; 1.7*theta2_i];
            theta_4_ = [theta_4_; theta4_i];
        end

%         dtheta_4 = cellfun(@(c) c{2}, results_data.vision_data);
%         ddtheta_4 = cellfun(@(c) c{3}, results_data.vision_data);
%         theta_4_ = cellfun(@(c) c{4}, results_data.vision_data);
%         theta_1_ = cellfun(@(c) c{5}, results_data.vision_data);
%         theta_2_ = cellfun(@(c) c{6}, results_data.vision_data);
        theta_1_mean = mean(theta_1_);
        theta_1_ = theta_1_ - theta_1_mean*ones(length(theta_1_),1);
        
        theta_2_mean = mean(theta_2_);
        theta_2_ = theta_2_ - theta_2_mean*ones(length(theta_2_),1);
        
        % Design a low-pass filter
        cutoffFreq = 8/20; % Cutoff frequency
        [b, a] = butter(3, cutoffFreq, 'low');
        %theta_4 = filtfilt(b, a, theta_4_);% - 0.0179*ones(length(theta_4),1);
        theta_1 = filtfilt(b, a, theta_1_);
        theta_2 = filtfilt(b, a, theta_2_);
        %%
        %theta_4 = imgaussfilt(theta_4_, 14);% - 0.0179*ones(length(theta_4),1);
        theta_4 = linear_kalman(ddtheta_4, [theta_4_, dtheta_4]);
        theta_4 = theta_4(:,1);

        ud = [ud; ddtheta_4];
        yd = [yd; theta_1, theta_2, theta_4, dtheta_4];
        wd = [wd; ud, yd];

        w = [acceleration, theta_1, theta_2, theta_4, velocity];
        H_w = blkhank(w, L);
        P = [P, H_w];

    end
end

%         theta_3 = pi/2 - 0.76;
%         rope_length = 1 - 0.2; 
%         theta_4_offset = 0;
%         z_boom = 1.57;
%         theta_4 = [];
%         theta_1 = [];
%         theta_2 = [];
%         for i=1:length(u)
%             [theta1_, theta2_, theta4_] = pixels_to_world(u(i), v(i), t4(i), theta_3, rope_length, theta_4_offset, z_boom);
%             theta_1 = [theta_1; theta1_];
%             theta_2 = [theta_2; theta2_];
%             theta_4 = [theta_4; theta4_];
%         end
%         theta_1_mean = mean(theta_1);
%         theta_1 = theta_1 - theta_1_mean*ones(length(theta_1),1);
%         % theta_2 = cellfun(@(c) c{6}, vision_data);
%         theta_2_mean = mean(theta_2);
%         theta_2 = theta_2 - theta_2_mean*ones(length(theta_2),1);
%         
%         % Design a low-pass filter
%         cutoffFreq = 0.06; % Cutoff frequency
%         [b, a] = butter(5, cutoffFreq, 'low');
%         theta_4 = filtfilt(b, a, theta_4);
%         theta_1 = filtfilt(b, a, theta_1);
%         theta_2 = filtfilt(b, a, theta_2);
