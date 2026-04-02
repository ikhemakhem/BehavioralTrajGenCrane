% Convert camera pixel measurements into crane load and boom angles.

function [theta_1, theta_2, theta_4] = pixels_to_world(u, v, theta_4, theta_3, rope_length, theta_4_offset, z_boom)
    %%% Given from Carne:
%     theta_3 = pi/2-deg2rad(6.2); %[rad]  vertical crane boom angle, always fixed.
    l = rope_length; %[m]  rope length always fixed.
    l_b = z_boom/cos(theta_3); %[m] boom length
    %theta_4 = deg2rad(theta_4);
    
    %%% Given from Camera Calibration:
    % Camera matrix K=[f_x, 0, c_x; 0, f_y, c_y; 0, 0, 1]
    K = [779.1104169208932, 0, 635.4402695731219;...
         0, 778.6639838714995, 374.9313955118369;...
         0, 0, 1];
    %%% Given from PnP Problem:
    % load all three SolvePnP results
%     load("Camera_Coordinates_1.mat","A_WP","W_r_OP");
%     R_0 = A_WP;
%     T_0 = W_r_OP;
%     load("Camera_Coordinates_2.mat","A_WP","W_r_OP");
%     R_0 = R_0 + A_WP;
%     T_0 = T_0 + W_r_OP;
%     load("Camera_Coordinates_3.mat","A_WP","W_r_OP");
%     R_0 = R_0 + A_WP;
%     T_0 = T_0 + W_r_OP;
%     % compute mean
%     A_WP = 1/3*R_0;     % rotation from camera (fixed in 0 degree) to world
%     W_r_OP = 1/3*T_0;   % coordinate of the camera (fixed in 0 degree) in world coordinates.

    load("Camera_Coordinates_June_Mean.mat","A_WP","W_r_OP");
    
    %%% theta_4 correction:
    theta_4 = theta_4_offset+ theta_4;
    %%% Height z
    %z = z_boom-l;

    z_camera = W_r_OP(3)-l;
    
    normalized_coord = (K\eye(3))*[u;v;1];  
    S_r_SI = z_camera*normalized_coord;  % S_r_SI: coordinates of point I expressed in the S (moving camera) coordinate system
    
    % Camera extrinsic parameters:
    A_SP = [cos(theta_4) -sin(theta_4) 0;
            sin(theta_4) cos(theta_4) 0;
            0 0 1];

    %%% World Coordinates:
    W_r_OI = A_WP*A_SP'*S_r_SI+A_SP*W_r_OP;
    
    %%% Minimal Coordinates:
    alpha = l_b*sin(theta_3)*cos(theta_4);
    beta = l_b*sin(theta_3)*sin(theta_4);
    load_angles = (1/l)* [cos(theta_4) sin(theta_4); -sin(theta_4) cos(theta_4)] * [W_r_OI(1)-alpha; W_r_OI(2)-beta];
    
    %%% Results
    theta_1 = load_angles(1);
    theta_2 = load_angles(2);

