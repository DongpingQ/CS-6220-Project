%% This algorithm finds SVD through SRFT by means of the interpolative decompostion
% Set-up

m = 30; % number of rows
n = 40; % number of columns
k = 5; % lower rank we wish to approximate
l = k+8; % common choice of l
A = rand(m,n); % A is the m*n target matrix we wish to find an approximation

%%
% Construct an interpolative decomposition of A, A = BP
% Use a rank k ID decomposition algorithm instead (see reference)
[~,Y] = srftmultback(l,A);
[Qy,Ry,Py]=qr(Y);
Qy1=Qy(:,1:k);
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