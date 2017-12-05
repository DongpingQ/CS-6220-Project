function [Q,X] = lowrankapprox(A,k,l,method,Gaussian_sketch,SRFT_sketch)
% This function performs a one-pass low-rank approximation to A with
% sketches of A
% Input: A: target matrix
%        k: a number slightly larger than r to perform the algorithm
%        l: a number slightly larger than k, details of k and l see referenced paper
%        method: a boolean value determining the random matrix used; 'G'
%                for Gaussian random matrix and 'S' for SRFT
%        Gaussian_sketch, SRFT_sketch, srftmultfront, srftmultback: helper functions
% Output: a low rank approximation of A, A = QX

if method == 'G'
    [Y,W,psi] = Gaussian_sketch(A,k,l);
else
    [Y,W,psi] = SRFT_sketch(A,k,l);
end

[Q,~]=qr(Y,0);
[U,T]=qr(psi*Q);
X=T\(U'*W);
end