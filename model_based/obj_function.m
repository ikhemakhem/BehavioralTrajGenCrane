function cost = obj_function(x,H_BoomPoints)
 
global Len_sv
 XH = 0.0;

 for e = 1:Len_sv
     XH       =  XH + (x(1,e)   -   H_BoomPoints(e))^2;
 end
  
    
WH        = 0.01;

WH        =1/(WH^2*Len_sv);


Dist  = WH*XH ;

cost = 10.0*x(1,3*Len_sv+1)  + 1*Dist;

% cost = x(1,3*Len_sv+1);