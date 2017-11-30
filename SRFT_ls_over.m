function x = SRFT_ls_over(A,b,eps,l,conj_grad,randH)
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
% T=srftmult_ls(l,H);

%%
% Least squares problem algorithm
% E=T*A;
E=S*fft(D*H*A);
[Q,R,Per]=qr(E,0);
iden=eye(n);
P=R*iden(Per,:);
% z=P\(Q'*(T*b));
z=P\(Q'*S*fft(D*H*b));
x=conj_grad(A,P,b,z,eps);
end