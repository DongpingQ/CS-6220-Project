function [C] = SRHT(A, s)
% This is the SRHT(Subsampled Randomized Hadarmard Transform) algorithm
% ======================================
% Input: A: matrix to be sketched, A is m by n
%        s: S = DHP/sqrt(sn), S is n by s
%           G is sampled i.i.d. from N(0,1)
% ======================================

n = size(A,2);
sgn = randi(2,[1,n])*2 - 3;    % P(1) = P(-1) = 0.5;  
A = bsxfun(@times, A, sgn);    % flip the sign of each column w.p. 0.5
n = 2^(ceil(log2(n)));
C = (fwht(A',n));              % fast Walsh - Hadarmard transform
idx = sort(randsample(n,s));
C = C(:,idx);
C = C*(n/sqrt(s));

end