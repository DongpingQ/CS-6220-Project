function [xtilde] = LSRN(A, b, C)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

Atilde = (A*C)';
[U, S, V] = svd(Atilde);
N = V\S;
y = pinv(A*N)*b;
xtilde = N*y;

end

