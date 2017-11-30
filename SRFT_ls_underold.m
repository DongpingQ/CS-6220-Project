function x = SRFT_ls_underold(A,b,eps,l,srftmult_ls,conj_grad,randH,SRFT_ls_overold)
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
T=srftmult_ls(l,H);

%%
% Least squares problem algorithm
S=T*A';
[Q,R]=qr(S);
t=R(1:m,1:m)'\b;
z=Q(:,1:m)*t;
c=T'*z;
y=SRFT_ls_overold(A',c,eps^2*l/(4*n),l,srftmult_ls,conj_grad,randH);
x=A'*y;
end