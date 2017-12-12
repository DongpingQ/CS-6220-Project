function [U,S,V] = SRFT_svd(A,k)
% This algorithm finds SVD through SRFT by means of the interpolative decompostion
% Input: A: target m-by-n matrix
%        k: the desired fixed low rank of A, k <= min(m,n)
%        srftmultback: helper function
% Output: The SVD of A, A=USV'
%%
% Set-up
[m,n]=size(A);
l = k+8; % common choice of l

%%
% Construct an interpolative decomposition of A, A = BP
% Use a rank k ID decomposition algorithm instead (see reference)
S=randi(m,1,l);
D=2*pi*rand(m,1);
D=cos(D)+1i*sin(D);
Y=fft(bsxfun(@times,A,D));
Y=Y(S,:);
% [Qy,Ry,Py]=qr(Y);
% Qy1=Qy(:,1:k);

[~,Ry,Py]=qr(Y,'vector');
iden=eye(n);
Perm=iden(Py,:);
Ry1=Ry(1:k,:);
Ry11=Ry1(:,1:k);
Ry12=Ry1(:,k+1:n);

B=A(:,Py(1:k));
P=Perm(1:k,:)+Ry11\(Ry12*Perm(k+1:n,:));

%%
% Construct the SVD of the matrix A = USV^{*}
[Q,Rp]=qr(P',0);
L=Rp';
C=B*L;
[U,S,W]=svds(C,k);
V=Q*W;
% Then A = USV^{*}