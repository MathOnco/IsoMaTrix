function [x_final, y_final] = UVW_to_XY(x)

L = 1;
H = (L/2)*tan(60*pi/180);

P = x(:,1);
Q = x(:,2);
x_final = (P./(tan(pi/3)) + (1-P-Q)./(sin(pi/3)))*sin(pi/3); 
y_final = P * H;
  
end
