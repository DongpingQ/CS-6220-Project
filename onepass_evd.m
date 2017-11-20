function [U,lam] = onepass_evd(omega,Y,Q)
% This function computes the eigenvalue decomposition of a Hermitian matrix
% A with a randomized algorithm without using A in this function
% Input: omega: random Gaussian matrix used in Stage A to find Q
%        Y: A*omega
%        Q: Orthonormal matrix the range finder algorithm finds A=AQQ'
% Output: U, lam: the eigenvalue decomposition of A, A=U*lam*U'

% Need a better least-squares solver
Bapprox=(Q'*Y)/(Q'*omega);

[V,lam]=eig(Bapprox);
U=Q*V;
end