function [Q,omega,Y] = subiter(A,q,l)
% This function performs a randomized subspace iteration algorithm to do a
% rnage finder to a matrix A with slow decaying spectrum
% Algorithm 4.4 in Halko et al. paper
% Input: A: target matrix
%        q: maximum step size
%        l: a number slightly larger than the desired rank k
% Output: Q: the approximate range matrix, A=AQQ'
%         omega: random Gaussian matrix used
%         Y: sample matrix Y=A*omega
[~,n]=size(A);
omega=normrnd(0,1,n,l);
Y0=A*omega;
[Q0,~]=qr(Y0);
for j=1:q
    Yhat=A'*Q0;
    [Qhat,~]=qr(Yhat);
    Y=A*Qhat;
    [Q1,~]=qr(Y);
    Q0=Q1;
end
Q=Q1;
end