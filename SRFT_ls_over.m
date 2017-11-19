function x = SRFT_ls_over(A,b,eps,k,srftmult,conj_grad,randH)
%% This algorithm solves overdetermined Least Squares Problem through SRFT
% Set-up

[m,n]=size(A);
l=k+8;

%%
H=randH(m);
T=srftmult(m,l,H);

%%
% Least squares problem algorithm
E=T*A;
[Q,R,Per]=qr(E,0);
iden=eye(n);
P=R*iden(:,Per);
z=P\(Q'*(T*b));
y=conj_grad(A,P,b,z,eps);
x=P\y;
end