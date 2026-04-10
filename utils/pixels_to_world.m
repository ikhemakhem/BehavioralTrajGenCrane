% Convert a camera pixel observation into crane load and boom angles.
%
% Uses the calibrated camera intrinsic matrix K and extrinsic parameters
% (rotation A_WP, translation W_r_OP) stored in data/Camera_Coordinates_June_Mean.mat.
% The calibration was obtained with a standard RGB camera using a PnP solver
% (see data/Camera_Coordinates_June_Mean.mat).
%
% Inputs:
%   u              - pixel column of the load marker
%   v              - pixel row of the load marker
%   theta_4        - current horizontal boom angle [rad]
%   theta_3        - fixed vertical boom elevation angle [rad]
%   rope_length    - rope length [m]
%   theta_4_offset - angular offset applied to theta_4 before back-projection [rad]
%   z_boom         - vertical height of the boom tip [m]
%
% Outputs:
%   theta_1 - load sway angle in the radial direction     [rad]
%   theta_2 - load sway angle in the tangential direction [rad]
%   theta_4 - horizontal boom angle (offset-corrected)    [rad]

function [theta_1, theta_2, theta_4] = pixels_to_world(u, v, theta_4, theta_3, rope_length, theta_4_offset, z_boom)
    % Load camera extrinsic calibration from the bundled data folder.
    thisFolder  = fileparts(mfilename('fullpath'));
    calibFile   = fullfile(fileparts(thisFolder), 'data', 'Camera_Coordinates_June_Mean.mat');
    calib       = load(calibFile, 'A_WP', 'W_r_OP');
    A_WP  = calib.A_WP;
    W_r_OP = calib.W_r_OP;

    l   = rope_length;
    l_b = z_boom / cos(theta_3);   % effective boom length

    % Apply theta_4 offset.
    theta_4 = theta_4_offset + theta_4;

    % Camera intrinsic matrix K (calibrated for the experimental setup).
    K = [779.1104169208932,  0,                  635.4402695731219;
         0,                  778.6639838714995,  374.9313955118369;
         0,                  0,                  1              ];

    % Back-project pixel to camera-frame ray, then scale to load height.
    z_camera        = W_r_OP(3) - l;
    normalized_coord = (K \ eye(3)) * [u; v; 1];
    S_r_SI          = z_camera * normalized_coord;

    % Camera extrinsic rotation (moving camera, rotates with boom).
    A_SP = [cos(theta_4), -sin(theta_4), 0;
            sin(theta_4),  cos(theta_4), 0;
            0,             0,            1];

    % World coordinates of the load.
    W_r_OI = A_WP * A_SP' * S_r_SI + A_SP * W_r_OP;

    % Project onto minimal crane coordinates.
    alpha       = l_b * sin(theta_3) * cos(theta_4);
    beta        = l_b * sin(theta_3) * sin(theta_4);
    load_angles = (1/l) * [ cos(theta_4),  sin(theta_4);
                            -sin(theta_4),  cos(theta_4)] * [W_r_OI(1) - alpha;
                                                              W_r_OI(2) - beta];
    theta_1 = load_angles(1);
    theta_2 = load_angles(2);
end
