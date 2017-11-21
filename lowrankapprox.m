function [Q,X] = lowrankapprox(A,k,l,method,Gaussian_sketch,SRFT_sketch)
% This function performs a low-rank approximation to A using only A once

if method == 'G'
    [Y,W,psi] = Gaussian_sketch(A,k,l);
else
    [Y,W,psi] = SRFT_sketch(A,k,l);
end

[Q,~]=qr(Y,0);
[U,T]=qr(psi*Q);
X=T\(U'*W);
end