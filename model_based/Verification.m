function x_dot=Verification(t,y,Vel,Acc)

    global m_load L g sv Len_sv th3 l

x_dot=zeros(4,1);

x_dot(1)=y(2);
x_dot(2)=( -l^2*(-2*y(4)*Vel - y(3)*Acc - y(1)*Vel^2 + y(1)*y(2)^2 + y(1)*y(4)^2) ...
          -L*l*(- Vel^2*sin(th3)) ...
          - g*l*y(1))/(l^2+l^2*y(1)^2);
x_dot(3)=y(4);

x_dot(4)= ( -L*l*(Acc*sin(th3) ) ...
            -l^2*(2*y(2)*Vel + y(1)*Acc + y(3)*y(4)^2 + y(2)^2*y(3) - y(3)*Vel^2) ...
            - g*l*y(3))/(l^2+l^2*y(3)^2);

x_dot=x_dot(:);
end



