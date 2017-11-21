function [Q,omega,Y] = gaussianrange_basic(A,l)
% This function performs a randomized fixed-rank range finder to a matrix A
% Algorithm 4.1 in Halko et al. paper
% Input: A: target matrix
%        l: a number slightly larger than the desired rank k
% Output: Q: the approximate range matrix, A=AQQ'
%         omega: random Gaussian matrix used
%         Y: sample matrix Y=A*omega

[~,n]=size(A);
omega=randn(n,l)+1i*randn(n,l);
Y=A*omega;
[Q,~]=qr(Y);
end