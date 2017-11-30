function [A] = kSVDMatrix(m, n, r, k)
% This is the k-SVD matrix generation algorithm
% ======================================
% Input: A: matrix to generate, A is m by n
%        k: first k singular value is dominant  
%        r: rank of A
% ======================================

s = zeros(r,1);
for i = 1:k
   s(i) = 10*(2*k+1-i)/k; 
end

for i = k+1:r
   s(i) = 5*(2*r+1-i)/(2*r-k); 
end

U = get_orthonormal(m,r);

V = get_orthonormal(n,r)';

A = zeros(m,n);
for i = 1:r
    A = A + s(i)*U(:,i).*V(i,:);
end

end

