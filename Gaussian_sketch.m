function [Y,W,psi] = Gaussian_sketch(A,k,l)
% This function performs a sketch of the target matrix A with Gaussian
% random matrices
% Consider complex cases since we need to compare it with SRFT
% Input: A: target matrix
%        k: a number slightly larger than r to perform the algorithm
%        l: a number slightly larger than k, details of k and l see referenced paper
% Output: Y: one sketch of A that Y = A*omega, omega is a random matrix
%         W: one sketch of A that W = psi*A, psi is another random matrix
%         psi: the random matrix used to compute W

[m,n]=size(A);
omega=randn(n,k)+1i*randn(n,k);
psi=randn(l,m)+1i*randn(l,m);
Y=A*omega;
W=psi*A;
end