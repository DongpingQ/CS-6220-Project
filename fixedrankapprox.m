function [Q,sigma,V] = fixedrankapprox(A,r,k,l,method,lowrankapprox,Gaussian_sketch,SRFT_sketch)
% This function performs SVD of A using A only once

[Q,X] = lowrankapprox(A,k,l,method,Gaussian_sketch,SRFT_sketch);
[U,sigma,V]=svds(X,r);
Q=Q*U;
end