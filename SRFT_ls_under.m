function x = SRFT_ls_under(A,b,eps,l,srftmult_ls,conj_grad,randH,SRFT_ls_over)
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
% Z1=2*pi*rand(n,1);
% Z2=2*pi*rand(n,1);
% Z1=cos(Z1)+1i*sin(Z1);
% Z2=cos(Z2)+1i*sin(Z2);
% Pi1=randperm(n);
% Pi2=randperm(n);
% 
% Y=bsxfun(@times,A',Z2);
% Y=Y(Pi2,:);
% tic
% for j=1:n-1
%     th=2*pi*rand;
%     givens=[cos(th) sin(th);-sin(th) cos(th)];
%     Y(n-j:n-j+1,:)=givens*Y(n-j:n-j+1,:);
%     z(n-j:n-j+1,:)=givens*z(n-j:n-j+1,:);
% end
% toc
% Y=bsxfun(@times,Y,Z1);
% Y=Y(Pi1,:);
% for j=1:n-1
%     th=2*pi*rand;
%     givens=[cos(th) sin(th);-sin(th) cos(th)];
%     Y((n-j):(n-j+1),:)=givens*Y((n-j):(n-j+1),:);
%     z((n-j):(n-j+1),:)=givens*z((n-j):(n-j+1),:);
% end
% E=fft(bsxfun(@times,Y,D));
Es=fft(bsxfun(@times,A',D));
Es=Es(S,:);
Es=Es/sqrt(m);
%T=srftmult_ls(l,H);

%%
% Least squares problem algorithm
%S=T*A';
[Q,R]=qr(Es);
t=R(1:m,1:m)'\b;
z=Q(:,1:m)*t;
%c=T'*z;
c=(S*fft(D*H))'*z;
y=SRFT_ls_over(A',c,eps^2*l/(4*n),l,srftmult_ls,conj_grad,randH);
x=A'*y;
end