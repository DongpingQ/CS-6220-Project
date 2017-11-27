m = 10;    n = 8;    r = 5;    k = 3;

s = zeros(r,1);
for i = 1:k
   s(i) = 10*(2*k+1-i)/k + rand; 
end

for i = k+1:r
   s(i) = 5*(2*k+1-i)/k + rand; 
end

U = get_orthonormal(m,r);

V = get_orthonormal(n,r)';

A = zeros(m,n);
for i = 1:r
    A = A + s(i)*U(:,i).*V(i,:);
end