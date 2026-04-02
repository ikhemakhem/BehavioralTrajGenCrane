function [c,ceq,pos_h,sp_h,ac_h,t_sw,y_sw]=nonlcon(x,Zero_Value,H_BoomPoints,...
                                          LB_Max_Sway,UB_Max_Sway,LB_MAX_Sway_Vel,UB_MAX_Sway_Vel, ...
                                          LB_Final_Sway,UB_Final_Sway,LB_Final_Vel_Sway,UB_Final_Vel_Sway)
    
  global m_load L g sv Len_sv th3 l


%% 
    Max_BH   = max(H_BoomPoints);
    Min_BH   = min(H_BoomPoints);


   Vel_Boom5_min = -0.1724;
   Vel_Boom5_max =  0.1724;

    Acc_Boom5_min =  -0.01724;
    Acc_Boom5_max =   0.01724;
    
%%
     ceq5_1  = zeros(1,Len_sv);  ceq5_2  = zeros(1,Len_sv);  ceq5_3  = zeros(1,Len_sv);
     ceq5_4  = zeros(1,Len_sv);  ceq5_5  = zeros(1,Len_sv);  ceq5_6  = zeros(1,Len_sv);

     c5_1  = zeros(1,Len_sv);   c5_2  = zeros(1,Len_sv);  c5_3  = zeros(1,Len_sv);
     c5_4  = zeros(1,Len_sv);   c5_5  = zeros(1,Len_sv);  c5_6  = zeros(1,Len_sv);
     
     
     c1_1  = zeros(1,Len_sv);      c2_1  = zeros(1,Len_sv);
     c1_2  = zeros(1,Len_sv);      c2_2  = zeros(1,Len_sv);
     c1_3  = zeros(1,Len_sv);      c2_3  = zeros(1,Len_sv);
     c1_4  = zeros(1,Len_sv);      c2_4  = zeros(1,Len_sv);
     

%%      
    ceq5_1(1,1)       = 1*(x(1,1)  - H_BoomPoints(1));
    ceq5_2(1,1)       = 1*x(1,Len_sv+1);
    ceq5_3(1,1)       = 1*x(1,2*Len_sv+1);

    for e = 2:Len_sv-1
        ceq5_2(1,e)   =  (((x(1,e)-x(1,e-1))/x(1,end))               -  x(1,Len_sv+e+1));
        ceq5_3(1,e)   =  (((x(1,Len_sv+e)-x(1,Len_sv+e-1))/x(1,end)) -  x(1,2*Len_sv+e+1));
    end
    
    ceq5_1(1,Len_sv)       = 10*(x(1,Len_sv)  - H_BoomPoints(end));
    ceq5_2(1,Len_sv)       = 1*x(1,2*Len_sv);
    ceq5_3(1,Len_sv)       = 1*x(1,3*Len_sv);



        c5_1(1,1:Len_sv)     = -min(x(1,1:Len_sv)) + Min_BH;
        c5_2(1,1:Len_sv)     =  max(x(1,1:Len_sv)) - Max_BH;

        c5_3(1,1:Len_sv)     = -min(x(1,Len_sv+1:2*Len_sv)) + Vel_Boom5_min;
        c5_4(1,1:Len_sv)     =  max(x(1,Len_sv+1:2*Len_sv)) - Vel_Boom5_max;

        c5_5(1,1:Len_sv)     = -min(x(1,2*Len_sv+1:3*Len_sv)) + Acc_Boom5_min;
        c5_6(1,1:Len_sv)     =  max(x(1,2*Len_sv+1:3*Len_sv)) - Acc_Boom5_max;


               
%%  %%   Solving ODE to satify the load-sway Constraints
ic =[0 0 0 0];
tt=sv*(Len_sv-1)*x(1,3*Len_sv+1);
pos_h = x(1,1:Len_sv);  
sp_h = x(1,Len_sv+1:2*Len_sv);
ac_h = x(1,2*Len_sv+1:3*Len_sv);

for i=1:Len_sv-1
        span=[tt(i) tt(i+1)];
        if tt(i)== tt(i+1)
            span=[tt(i) tt(i+1)+0.01];
        end
        [t,y] = ode45(@(t,y)Verification_Model(t,y, sp_h(i),ac_h(i)), span, ic);
        ic=y(end,:);
        t_sw(i+1,1)=t(end);
        y_sw(i+1,:)=y(end,:);
    end

     c1_1(1,1:Len_sv-1)     =  -min(y_sw(:,1)) + LB_Max_Sway;
     c1_2(1,1:Len_sv-1)     =  max(y_sw(:,1)) - UB_Max_Sway;
     c1_3(1,1:Len_sv-1)     = -min(y_sw(:,2)) + LB_MAX_Sway_Vel;
     c1_4(1,1:Len_sv-1)     =  max(y_sw(:,2)) - UB_MAX_Sway_Vel;

     c2_1(1,1:Len_sv-1)     = -min(y_sw(:,3)) + LB_Max_Sway;
     c2_2(1,1:Len_sv-1)     =  max(y_sw(:,3)) - UB_Max_Sway;
     c2_3(1,1:Len_sv-1)     = -min(y_sw(:,4)) + LB_MAX_Sway_Vel;
     c2_4(1,1:Len_sv-1)     =  max(y_sw(:,4)) - UB_MAX_Sway_Vel;

     c1_1(1,Len_sv)         =  1*(-y_sw(end,1) + LB_Final_Sway);
     c1_2(1,Len_sv)         =  1*(y_sw(end,1) - UB_Final_Sway);
     c1_3(1,Len_sv)         =  1*(-y_sw(end,2) + LB_Final_Vel_Sway);
     c1_4(1,Len_sv)         =  1*(y_sw(end,2) - UB_Final_Vel_Sway);
     
     c2_1(1,Len_sv)         =  1*(-y_sw(end,3) + LB_Final_Sway);
     c2_2(1,Len_sv)         =  1*(y_sw(end,3) - UB_Final_Sway);
     c2_3(1,Len_sv)         =  1*(-y_sw(end,4) + LB_Final_Vel_Sway);
     c2_4(1,Len_sv)         =  1*(y_sw(end,4) - UB_Final_Vel_Sway);

    
  %%

    ceq = [ceq5_1 ceq5_2 ceq5_3];
    c   = [c5_3 c5_4  c5_5 c5_6 ...
           c1_1  c1_2 c1_3 c1_4 ...
           c2_1  c2_2  c2_3  c2_4];





