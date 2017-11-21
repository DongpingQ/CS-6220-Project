clear all 
format compact

m = 10;    n = 100;    eps = 1.0e-1;

s = floor(m/eps^2);

A = 10*rand(m, n) + 10*sin(rand(m, n));

C = GaussianProjection(A, s);

Cinv = pinv(C, eps);

X = A - C*Cinv*A;

nCfro = norm(X,'fro');
fprintf('Frobenius norm: ')
disp(nCfro)