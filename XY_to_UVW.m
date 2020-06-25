%% reads in [x, y] pair to return [x1 x2 x3] coords
function [P,Q,w] = XY_to_UVW(p)
    L = 1;
    H = (L/2)*tan(60*pi/180);
   
    y1 = p(:,1);
    y2 = p(:,2);
    P = y2 ./ H;

    temp = (P-1)*sin(pi/3) + (y1./(sin(pi/3))) - P*tan(pi/3);
    Q = -temp ./ sin(pi/3);
    
    w = 1-P-Q;
    
end