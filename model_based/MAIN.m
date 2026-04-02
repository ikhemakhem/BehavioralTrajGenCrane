    clear all
    close all
    clc


    tic
%%
    global m_load L g sv Len_sv th3 l
    g=9.81;
    m_load=3; 
    L=2.2511;

% Original

    th3 =  0.8108;
%     H_BoomPoints = [0 3 6 9 12 15 18 21 24 27 30 33 36 39 42 45]*pi/180;
    H_BoomPoints = [45 48 51 54 57 60 63 66 69 72 75 78 81 84 87 90]*pi/180;
    l   =  1.0;
    
%%
    Zero_Value    = 0;
    
    Len_sv        = length(H_BoomPoints);
    sv            = linspace(0,1,Len_sv);
    
  
%%
% Swaying Angles Condition
 
    Max_Sway           = 2.0*pi/180;
    MAX_Sway_Vel       = 1.0*pi/180;
    LB_Max_Sway        = -Max_Sway; 
    UB_Max_Sway        =  Max_Sway;
    LB_MAX_Sway_Vel    = -MAX_Sway_Vel;  
    UB_MAX_Sway_Vel    =  MAX_Sway_Vel;
         
    Final_Sway         = 0.1*pi/180;
    Final_Vel_Sway     = 0.5*pi/180;
    LB_Final_Sway      = -Final_Sway;
    UB_Final_Sway      =  Final_Sway;
    LB_Final_Vel_Sway  = -Final_Vel_Sway;
    UB_Final_Vel_Sway  =  Final_Vel_Sway;
    
    %%

    LB                              = zeros(1,3*Len_sv+1);
    LB(1,3*Len_sv+1)                = 0;
    LB(1,1:Len_sv)                  = -inf;
    LB(1,Len_sv+1:2*Len_sv)         = -0.1724;
    LB(1,2*Len_sv+1:3*Len_sv)       = -0.01724;
    
     
    UB                              = zeros(1,3*Len_sv+1);
    UB(1,3*Len_sv+1)                = inf;
    UB(1,1:Len_sv)                  = inf;
    UB(1,Len_sv+1:2*Len_sv)         = 0.1724;
    UB(1,2*Len_sv+1:3*Len_sv)       = 0.01724;


     
    x0                              = zeros(1,3*Len_sv+1);
    x0(1,3*Len_sv+1)                = 1.5;
    x0(1,1:Len_sv)                  = [45 48 51 54 57 60 63 66 69 72 75 78 81 84 87 90]*pi/180;
    x0(1,Len_sv+1:2*Len_sv)         = 0.1724;
    x0(1,2*Len_sv+1:3*Len_sv)       = 0.01724;



%%

    A=[];
    B=[];
    Aeq=[];
    Beq=[];
       
    options=optimoptions('fmincon', 'Algorithm', 'sqp', 'Display',  'iter-detailed', ...
        'MaxFunctionEvaluation', 100000, 'MaxIterations', 1000, ...
        'FunctionTolerance', 1e-10);



    [x,fval,exitflag,output]=fmincon(@(x)obj_function(x,H_BoomPoints), x0, A, B, Aeq, Beq, LB, UB, ... 
                                          @(x) nonlcon(x,Zero_Value,H_BoomPoints,...
                                          LB_Max_Sway,UB_Max_Sway,LB_MAX_Sway_Vel,UB_MAX_Sway_Vel, ...
                                          LB_Final_Sway,UB_Final_Sway,LB_Final_Vel_Sway,UB_Final_Vel_Sway), options);
                                             
    toc
    [c,ceq,pos_h,sp_h,ac_h,t_sw,y_sw]=nonlcon(x,Zero_Value,H_BoomPoints,...
                                          LB_Max_Sway,UB_Max_Sway,LB_MAX_Sway_Vel,UB_MAX_Sway_Vel, ...
                                          LB_Final_Sway,UB_Final_Sway,LB_Final_Vel_Sway,UB_Final_Vel_Sway);
    %%
 
    Time  =sv*(Len_sv-1)*x(1,end);
    pos_h = x(1,1:Len_sv);
    sp_h  = x(1,Len_sv+1:2*Len_sv);
    ac_h  = x(1,2*Len_sv+1:3*Len_sv);


    %%
    fig('fontsize',14,'width',18,'height',16,'units','centimeters')
    figure(1)
    subplot(3,1,1)
    plot(Time,pos_h*180/pi,'blue','LineWidth',2)
    ylabel({'$\theta_4$';'[deg.]'},'interpreter','latex')

    subplot(3,1,2)
    plot(Time,sp_h*180/pi,'blue','LineWidth',2)
    ylabel({'$\dot\theta_4$';'[deg./s]'},'interpreter','latex')

    subplot(3,1,3)
    plot(Time,ac_h*180/pi,'blue','LineWidth',2)
    ylabel({'$\ddot\theta_4$';'[deg./s]'},'interpreter','latex')


%%
fig('fontsize',14,'width',18,'height',16,'units','centimeters')
    figure(2)
    subplot(2,2,1)
    plot(Time,y_sw(:,1)*180/pi,'blue','LineWidth',2)
    ylabel({'$\theta_1$';'[deg.]'},'interpreter','latex')
  
     subplot(2,2,3)
    plot(Time,y_sw(:,2)*180/pi,':green','LineWidth',2)
    ylabel({'$\dot\theta_1$';'[deg./sec]'},'interpreter','latex')

    subplot(2,2,2)
    plot(Time,y_sw(:,3)*180/pi,'blue','LineWidth',2)
    ylabel({'$\theta_2$';'[deg.]'},'interpreter','latex')

    subplot(2,2,4)
    plot(Time,y_sw(:,4)*180/pi,':green','LineWidth',2)
    ylabel({'$\dot\theta_2$';'[deg./sec]'},'interpreter','latex')
    
    %%    
    col_header={'TIME','THETA_4','VEL_4','Acc_4','THETA_1','VEL_1','THETA_2','VEL_2'};
    xlswrite('Result.xlsx',col_header,'Sheet1','A1');
    xlswrite('Result.xlsx',[Time(:),pos_h(:),sp_h(:),ac_h(:),y_sw(:,1),y_sw(:,2),y_sw(:,3),y_sw(:,4)],'Sheet1','A2');

    x
    exitflag

    
