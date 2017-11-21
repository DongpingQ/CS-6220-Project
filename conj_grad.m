function soln = conj_grad(A,S,b,x,eps)
% This function computes the preconditioned conjugate gradient iterations
% given general A, b, precondition matrix S, initial guess x and tolerance
% eps, not a lot of iterations will be needed if the condition number of
% the preconditioned system is small
% Input: A: target matrix
%        S: precondition matrix
%        b: result vector
%        x: initial guess vector
%        eps: tolerance
% Output: soln: final solution after iterations
r=b-A*x;
p=S'\(A'*r);
s=p;
gam=norm(s)^2;
while gam>eps
    t=S\p;
    q=A*t;
    alpha=gam/(norm(q)^2);
    x=x+alpha*t;
    r=r-alpha*q;
    s=S'\(A'*r);
    gamnew=norm(s)^2;
    beta=gamnew/gam;
    gam=gamnew;
    p=s+beta*p;
end
soln=x;
end