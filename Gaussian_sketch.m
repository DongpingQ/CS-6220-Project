function [Y,W,psi] = Gaussian_sketch(A,k,l)
% This function performs a sketch of the target matrix A
% Consider complex cases since we need to compare it with SRFT

omega=randn(n,k)+1i*randn(n,k);
psi=randn(l,m)+1i*randn(l,m);
Y=A*omega;
W=psi*A;
end