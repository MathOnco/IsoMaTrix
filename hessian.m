function[D,lamda,V]=hessian(x,A)
%x is a 1*3 vector of the fixed point, B is the payoff matrix
%D is the 2*2 hessian matrix if we only use x1 and x2 as variable.
%lamba is the 2*2 matrix with diagonal element as the eigen values.
%V is a 3*2 matrix and each column is the eigen vector.

D=zeros(2,2);
V=zeros(3,2);


for i=1:2
    for j=1:2
        if i==j
            D(i,j)=A(i,:)*x'-x*A*x'+(A(i,i)-A(i,:)*x'-x*A(:,i))*x(i);
            D(i,j)=D(i,j)-(A(i,3)-A(3,:)*x'-x*A(:,3))*x(i);
        end
        
        if i~=j
            D(i,j)=(A(i,j)-A(j,:)*x'-x*A(:,j))*x(i);
            D(i,j)=D(i,j)-(A(i,3)-A(3,:)*x'-x*A(:,3))*x(i);
        end
    end
end

[v,lamda]=eig(D); 
V(1:2,1:2)=v;

V(3,1)=0-v(2,1)-v(1,1);  % the original eigen vector is is 2*1, we add extra element to make it 3*1 and sum up to 0.
V(3,2)=0-v(2,2)-v(1,2);



end