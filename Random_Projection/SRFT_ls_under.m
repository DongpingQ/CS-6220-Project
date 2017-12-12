function x = SRFT_ls_under(A,b,l,conj_grad,SRFT_ls_over)
% This algorithm solves underdetermined Least Squares Problem Ax=b with
% A of size m by n and m < n by means of SRFT
% Input: A: target matrix
%        b: target result vector
%        eps: tolerance for conjugate gradient iterations
%        l: parameter used in the algorithm n >= l > m
%        srftmultback, conj_grad, randH: helper functions
% Output: x: the solution to the overdetermined system
%%
% Set-up
[m,n]=size(A);

%%
% Construct another random matrix H to use together with SRFT
% H = th1*perm1*z1*th2*perm2*z2 and compute T = SRFT*H
S=randi(n,1,l);
D=2*pi*rand(n,1);
D=cos(D)+1i*sin(D);

Sh=zeros(n,l);
for j=1:l
    Sh(S(j),j)=1;
end

% Z1=2*pi*rand(n,1);
% Z2=2*pi*rand(n,1);
% Z1=cos(Z1)+1i*sin(Z1);
% Z2=cos(Z2)+1i*sin(Z2);
% Pi1=randperm(n);
% Pi2=randperm(n);
% [~,iPi1]=sort(Pi1);
% [~,iPi2]=sort(Pi2);
% 
% Y=bsxfun(@times,A',Z2);
% Y=Y(Pi2,:);
% 
% th=zeros(1,n-1);
% for j=1:n-1
%     th(j)=2*pi*rand;
%     givens=[cos(th(j)) sin(th(j));-sin(th(j)) cos(th(j))];
%     Y(n-j:n-j+1,:)=givens*Y(n-j:n-j+1,:);
% end
% 
% Y=bsxfun(@times,Y,Z1);
% Y=Y(Pi1,:);
% tht=zeros(1,n-1);
% for j=1:n-1
%     tht(j)=2*pi*rand;
%     givens=[cos(tht(j)) sin(tht(j));-sin(tht(j)) cos(tht(j))];
%     Y((n-j):(n-j+1),:)=givens*Y((n-j):(n-j+1),:);
% end
% Es=fft(bsxfun(@times,Y,D));

Es=fft(bsxfun(@times,A',D));
Es=Es(S,:);
Es=Es/sqrt(n);



%%
% Least squares problem algorithm
[Q,R]=qr(Es);
t=R(1:m,1:m)'\b;
z=Q(:,1:m)*t;

c=conj(D).*ifft(Sh*z)*n;
c=c/sqrt(n);

% for j=1:n-1
%     givens=[cos(tht(n-j)) -sin(tht(n-j));sin(tht(n-j)) cos(tht(n-j))];
%     c(j:j+1,:)=givens*c(j:j+1,:);
% end
% c=c(iPi1,:);
% c=conj(Z1).*c;
% for j=1:n-1
%     givens=[cos(th(n-j)) -sin(th(n-j));sin(th(n-j)) cos(th(n-j))];
%     c(j:j+1,:)=givens*c(j:j+1,:);
% end
% c=c(iPi2,:);
% c=conj(Z2).*c;

[y,~,~]=SRFT_ls_over(A',c,(1e-14*cond(A))^2*m/n,l,conj_grad);
x=A'*y;
end