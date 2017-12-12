function Q = Gaussianrange(A,eps)
% This function performs an adaptive randomized fixed-precision range finder to A
% through a fixed-rank algorithm
% Algorithm 4.2 in Halko et al. paper
% Input: A: target matrix
%        eps: precision desired that ||A-AQQ'|| \le eps
% Output: Q: range matrix with property specified above
[m,n]=size(A);
r=10;
w=randn(n,r)+1i*randn(n,r);
y=A*w;
j=0;
Qold=zeros(m,0);
while max(vecnorm(y(:,j+1:j+r))) > eps/(10*sqrt(2/pi))
    j=j+1;
    y(:,j)=(eye(m)-Qold*Qold')*y(:,j);
    q=y(:,j)/norm(y(:,j));
    Qnew=[Qold,q];
    wnew=randn(n,1)+1i*randn(n,1);
    ynew=(eye(m)-Qnew*Qnew')*A*wnew;
    y=[y,ynew];
    for i=j+1:j+r-1
        y(:,i)=y(:,i)-q*(q'*y(:,i));
    end
    Qold=Qnew;
end
Q=Qnew;
end