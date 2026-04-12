% Build the data matrix from the bundled processed crane-data files.
%
% Loads all paired input/result sequences from data/training_data/ and assembles
% the block Hankel matrix P together with aggregate data arrays ud, yd, wd.
% The data folder is located relative to this file, so the function works
% regardless of the caller's working directory.
%
% Input:
%   L  - trajectory window length (number of time steps per Hankel column)
%
% Outputs:
%   P  - (q*L x N_total) block Hankel data matrix
%   ud - (T_total x 1)   concatenated input signals   (boom acceleration ddtheta4)
%   yd - (T_total x 4)   concatenated output signals  [theta1, theta2, theta4, dtheta4]
%   wd - (T_total x 5)   concatenated stacked signals [ud, yd]

function [P, ud, yd, wd] = pageMatrixFromData(L)
    thisFolder = fileparts(mfilename('fullpath'));
    dataFolder = fullfile(fileparts(thisFolder), 'data', 'training_data');

    inputFiles = dir(fullfile(dataFolder, 'input_seq_*.mat'));
    inputFiles = sort({inputFiles.name});   % ensure sorted order

    P  = [];
    ud = [];
    yd = [];
    wd = [];

    for k = 1:length(inputFiles)
        % Load input sequence (variable: acceleration)
        inputData = load(fullfile(dataFolder, inputFiles{k}));
        ud_i = inputData.acceleration;

        % Load corresponding results sequence (variable: yd_i)
        resultsFile = sprintf('results_seq_%02d.mat', k);
        resultsData = load(fullfile(dataFolder, resultsFile));
        yd_i = resultsData.yd_i;

        ud   = [ud;   ud_i];         %#ok<AGROW>
        yd   = [yd;   yd_i];         %#ok<AGROW>
        wd_i = [ud_i, yd_i];
        wd   = [wd;   wd_i];         %#ok<AGROW>

        P = [P, blkhank(wd_i, L)];   %#ok<AGROW>
    end
end
