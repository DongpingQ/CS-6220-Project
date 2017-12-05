function [Y,W,psi] = SRFT_sketch(A,k,l)
% This function performs a sketch of A with SRFT
% Input: A: target matrix
%        k: a number slightly larger than r to perform the algorithm
%        l: a number slightly larger than k, details of k and l see referenced paper
%        srftmultfront, srftmultback: helper functions
% Output: Y: one sketch of A that Y = A*omega, omega is a random matrix
%         W: one sketch of A that W = psi*A, psi is another random matrix
%         psi: the random matrix used to compute W

[m,n]=size(A);
S2=randi(m,1,l);

% D2=2*pi*rand(m,1);
% D2=cos(D2)+1i*sin(D2);
% W=fft(bsxfun(@times,A,D2))/sqrt(m);
% W=W(S2,:);
% psi=fft(diag(D2))/sqrt(m);
% psi=psi(S2,:);
r=rand(1,m);
W=A;
W(r<1/2,:)=-W(r<1/2,:);
W=fft(W)/sqrt(m);
W=W(S2,:);
temp=ones(1,m);
temp(r<1/2)=-1;
psi=fft(diag(temp))/sqrt(m);
psi=psi(S2,:);

S1=randi(n,1,k);

% D1=2*pi*rand(1,n);
% D1=cos(D1)+1i*sin(D1);
% Y=fft(bsxfun(@times,A,D1))/sqrt(n);
% Y=Y(:,S1);
s=rand(1,n);
Y=A;
Y(:,s<1/2)=-Y(:,s<1/2);
Y=fft(Y)/sqrt(n);
Y=Y(:,S1);
end