function x = SRFT_ls_under(A,b,eps,k,srftmult,conj_grad,randH)
%% This algorithm solves underdetermined Least Squares Problem through SRFT
% Set-up

si=size(A);
m=si(1);
n=si(2);
l=k+8;

%%
% Construct another random matrix H to use together with SRFT
% H = th1*perm1*z1*th2*perm2*z2 and compute T = SRFT*H
H=randH(n);
T=srftmult(n,l,H);

%%
% Least squares problem algorithm
S=T*A';
[Q,R]=qr(S);
t=R(1:m,1:m)'\b;
z=Q(:,1:m)*t;
c=T'*z;
y=SRFT_ls_over(A',c,eps,k,srftmult,conj_grad);
x=A'*y;
end