    clear all
    close all
    clc
    
    
    data_csv_dtheta4=csvread('dtheta4.csv');
    M=transpose(data_csv_dtheta4);
    save('dtheta4.mat','M');


%     data_csv_velocity=csvread('velocity.csv');
%     M=transpose(data_csv_velocity);
%     save('velocity.mat','M');