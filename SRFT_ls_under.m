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
H=randH(n);
S = zeros(l,m);
p = randperm(m,l);
for j = 1:l
    S(j,p(j))=1;
end

% Generate D
d=zeros(1,m);
for j = 1:m
    theta=2*pi*rand;
    d(j)=cos(theta)+sin(theta)*1i;
end
D=diag(d);
%T=srftmult_ls(l,H);

%%
% Least squares problem algorithm
%S=T*A';
Z=S*fft(D*H*A');
[Q,R]=qr(Z);
t=R(1:m,1:m)'\b;
z=Q(:,1:m)*t;
%c=T'*z;
c=(S*fft(D*H))'*z;
y=SRFT_ls_over(A',c,eps^2*l/(4*n),l,srftmult_ls,conj_grad,randH);
x=A'*y;
end