function [x] = LSRN(A, b, C, tol)
% This is the LSRN algorithm
% ======================================
% Input: A: matrix to generate, A is m by n
%        k: first k singular value is dominant  
%        r: rank of A
% ======================================

Y = (C'*A)';
[Q, R] = qr(Y,0);
% T = inv(R);
Atilde = A/R;

x = lsqr(Atilde, b, tol);

end

