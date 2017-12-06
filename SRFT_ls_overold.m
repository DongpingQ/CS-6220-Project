function x = SRFT_ls_overold(A,b,eps,l,srftmultback,conj_grad,randH)
% This algorithm solves overdetermined Least Squares Problem Ax=b with
% A of size m by n and m >= n by means of SRFT
% Input: A: target matrix
%        b: target result vector
%        eps: tolerance for conjugate gradient iterations
%        l: parameter used in the algorithm m >= l > n
%        srftmultback, conj_grad, randH: helper functions
% Output: x: the solution to the overdetermined system

%%
% Set-up and construct random matrix
[m,n]=size(A);

H=randH(m);
[~,T]=srftmultback(l,H);

%%
% Least squares problem algorithm
E=T*A;
[Q,R,Per]=qr(E,0);
iden=eye(n);
P=R*iden(Per,:);
z=P\(Q'*(T*b));
x=conj_grad(A,P,b,z,eps);
end