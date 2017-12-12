function [U, Sigma, V] = OnePassSVD(A, k, Omega, Psi)
%This is the method of Using 

Y = A*Omega;
W = Psi*A;

[Q,~] = qr(Y,0);
[Q1,R1] = qr(Psi*Q);
X = R1\(Q1'*W);

[U,Sigma,V] = svds(X,k);

U = Q*U;
end

