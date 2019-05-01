function [y,t1,t2,d,t_vmax,v_max] = h_s(x,theta1,theta2,theta3,dm)
t1 = theta2*(-log(0.5+0.5*dm))^(1/theta3);
t2 = theta2*(-log(0.5-0.5*dm))^(1/theta3);
d = t2-t1;
X = x+t1;
X(X<0) = 0;
y = theta1*(1-exp(-(X/theta2).^theta3));
t_vmax = theta2*((theta3-1)/theta3)^(1/theta3)-t1;
v_max = theta1*theta3/theta2*exp((1-theta3)/theta3)*((theta3-1)/theta3)^((theta3-1)/theta3);
end

