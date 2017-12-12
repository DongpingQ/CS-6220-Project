function M = srftmultfront(k,A)
% This function generates a n*l sized subsampled random Fourier transform
% R=DFS and computes the matrix multiplication AR, R is formatted this way
% so that the matrix A is multiplied in the front
% Input: k: a parameter specifying the column number of the SRFT matrix
%        A: target matrix
% Output: M: the product of A and R, M = A*R

[~,n]=size(A);
% Generate S
S = zeros(n,k);
for j = 1:k
    sj=randi(n);
    S(sj,j)=1;
end

% Generate D
d=zeros(1,n);
for j = 1:n
    theta=2*pi*rand;
    d(j)=cos(theta)+sin(theta)*1i;
end
D=diag(d);

M = fft(A*D)*S;
end