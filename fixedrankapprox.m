function [Q,sigma,V] = fixedrankapprox(A,r,k,l,method,lowrankapprox,Gaussian_sketch,SRFT_sketch,srftmultfront,srftmultback)
% This function performs a one-pass SVD to matrix A
% Input: A: target matrix
%        r: specified lower rank r <= min(m,n)
%        k: a number slightly larger than r to perform the algorithm
%        l: a number slightly larger than k, details of k and l see referenced paper
%        method: a boolean value determining the random matrix used; 'G'
%                for Gaussian random matrix and 'S' for SRFT
%        lowrankapprox, Gaussian_sketch, SRFT_sketch, srftmultfront, srftmultback: helper functions
% Output: we approximate the SVD of A with Q*sigma*V'

[Q,X] = lowrankapprox(A,k,l,method,Gaussian_sketch,SRFT_sketch,srftmultfront,srftmultback);
[U,sigma,V]=svds(X,r);
Q=Q*U;
end