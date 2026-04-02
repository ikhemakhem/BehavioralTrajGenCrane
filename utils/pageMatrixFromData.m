% Build a data matrix from the bundled processed crane-data files.
function [P, ud, yd, wd] = pageMatrixFromData(L)
    % PAGE MATRIX FROM DATA
    % This function loads pre-saved input and result files numbered sequentially
    % and constructs the block Hankel matrix P and aggregated data matrices ud, yd, wd.

    % Move to the folder containing the data files
    thisFile  = mfilename('fullpath');        % full path of this .m file
    thisFolder = fileparts(thisFile);         % folder where this file is located
    parentFolder = fileparts(thisFolder);     % go one folder up
    dataFolder = fullfile(parentFolder, 'data', 'crane_data');  % bundled processed crane-data folder

    cd(dataFolder);

    % Initialize outputs
    P = [];
    ud = [];
    yd = [];
    wd = [];

    % Find all input files with the pattern input_seq_XX.mat
    inputFiles = dir('input_seq_*.mat');
    inputFiles = sort({inputFiles.name});  % ensure sorted order

    % Loop over each input file
    for k = 1:length(inputFiles)
        % Load input sequence
        inputData = load(inputFiles{k});   % expects variable 'ud_i'
        ud_i = inputData.ud;

        % Load corresponding results sequence
        resultsFile = sprintf('results_seq_%02d.mat', k); % match numbering
        resultsData = load(resultsFile);   % expects variable 'yd_i'
        yd_i = resultsData.yd;

        % Append to aggregate matrices
        ud = [ud; ud_i];
        yd = [yd; yd_i];
        wd_i = [ud_i, yd_i];
        wd = [wd; wd_i];

        % Construct block Hankel matrix
        P = [P, blkhank(wd_i, L)];
    end
end


