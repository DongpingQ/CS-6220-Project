function [U, S, V] = OnePassEVD(A, Omega)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

Y = A*Omega;
[Q,R] = qr(Y);
% using least square to find B_approx

[V,S] = eig(B_approx);
U = Q*V;

end

