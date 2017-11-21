function [C] = GaussianProjection(A, s)
% This is the Gaussian Projection algorithm
% ======================================
% Input: A: matrix to be sketched, A is m by n
%        s: S = G/sqrt(s), S is n by s
%           G is sampled i.i.d. from N(0,1)
% ======================================

n = size(A,2);
S = randn(n, s)/sqrt(s);
C = A*S;

end