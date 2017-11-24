function [Utilde, Stilde, Vtilde] = kSVDPrototype(A, k, C)
% This is the k-SVD prototype algorithm
% ======================================
% Input: A: matrix to do k-SVD, A is m by n
%        k: k-SVD parameter
%        C: sketch matrix: Gaussian or CountSketch 
% ======================================

[Q,R] = qr(C, 0);
[Ubar, Stilde, Vtilde] = svds(Q'*A, k);
Utilde = Q*Ubar;

end

