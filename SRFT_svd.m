function [U,S,V] = SRFT_svd(A,k,srftmultback)
% This algorithm finds SVD through SRFT by means of the interpolative decompostion
% Input: A: target m-by-n matrix
%        k: the desired fixed low rank of A, k <= min(m,n)
%        srftmultback: helper function
% Output: The SVD of A, A=USV'
%%
% Set-up
[~,n]=size(A);
l = k+8; % common choice of l

%%
% Construct an interpolative decomposition of A, A = BP
% Use a rank k ID decomposition algorithm instead (see reference)
[~,Y] = srftmultback(l,A);
% [Qy,Ry,Py]=qr(Y);
% Qy1=Qy(:,1:k);
[~,Ry,Py]=qr(Y);
Ry1=Ry(1:k,:);
Ry11=Ry1(:,1:k);
Ry12=Ry1(:,k+1:n);

P=[eye(k) Ry11\Ry12]*Py';
B=A*Py(:,1:k);

%%
% Construct the SVD of the matrix A = USV^{*}
[Q,Rp]=qr(P');
L=Rp';
C=B*L;
[U,S,W]=svd(C);
V=Q*W;
% Then A = USV^{*}