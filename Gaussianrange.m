function Q = Gaussianrange(A,eps)
[m,n]=size(A);
r=10;
w=normrnd(0,1,n,r);
y=A*w;
j=0;
Qold=zeros(m,0);
while max(vecnorm(y(:,j+1:j+r))) > eps/(10*sqrt(2/pi))
    j=j+1;
    y(:,j)=(eye(m)-Qold*Qold')*y(:,j);
    q=y(:,j)/norm(y(:,j));
    Qnew=[Qold,q];
    wnew=normrnd(0,1,n,1);
    ynew=(eye(m)-Qnew*Qnew')*A*wnew;
    y=[y,ynew];
    for i=j+1:j+r-1
        y(:,i)=y(:,i)-q*(q'*y(:,i));
    end
    Qold=Qnew;
end
Q=Qnew;
end