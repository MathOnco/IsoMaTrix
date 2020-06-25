function xdot = replicator(t,x,A)    
    f1 = A(1,1)*x(1) + A(1,2)*x(2) + A(1,3)*x(3);
    f2 = A(2,1)*x(1) + A(2,2)*x(2) + A(2,3)*x(3);
    f3 = A(3,1)*x(1) + A(3,2)*x(2) + A(3,3)*x(3);
    phi = f1.*x(1) + f2.*x(2) + f3.*x(3);
    
    xdot(1,1) = x(1)* ( f1 - phi ) ; 
    xdot(2,1) = x(2)* ( f2 - phi ) ;
    xdot(3,1) = x(3)* ( f3 - phi ) ;
    
end