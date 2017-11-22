function [C] = CountSketch(A, s)
% This is the Count SKetch algorithm
% ======================================
% Input: A: matrix to be sketched, A is m by n
%        s: S = sparse(1,-1,p(1) = p(-1) = 0.5), S is n by s
%           G is sampled i.i.d. from N(0,1)
% ======================================

[m, n] = size(A);
sgn = randi(2,[1,n])*2 - 3;    % P(1) = P(-1) = 0.5;  
A = bsxfun(@times, A, sgn);    % flip the sign of each column w.p. 0.5
ll = randsample(s, n, true);   % sample n items from [s] with replacement
C = zeros(m,s);

for j = 1:n
    C(:,ll(j)) = C(:,ll(j)) + A(:,j);
end

end