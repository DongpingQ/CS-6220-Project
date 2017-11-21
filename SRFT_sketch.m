function [Y,W,psi] = SRFT_sketch(A,k,l,srftmultfront,srftmultback)
% This function performs a sketch of A with SRFT
% Input: A: target matrix
%        k: a number slightly larger than r to perform the algorithm
%        l: a number slightly larger than k, details of k and l see referenced paper
%        srftmultfront, srftmultback: helper functions
% Output: Y: one sketch of A that Y = A*omega, omega is a random matrix
%         W: one sketch of A that W = psi*A, psi is another random matrix
%         psi: the random matrix used to compute W

Y=srftmultfront(k,A);
[psi,W]=srftmultback(l,A);
end