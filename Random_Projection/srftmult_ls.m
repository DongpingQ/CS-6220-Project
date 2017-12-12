function M = srftmult_ls(l,A)
[m,~]=size(A);
% Generate S
S = zeros(l,m);
p = randperm(m,l);
for j = 1:l
    S(j,p(j))=1;
end

% Generate D
d=zeros(1,m);
for j = 1:m
    theta=2*pi*rand;
    d(j)=cos(theta)+sin(theta)*1i;
end
D=diag(d);

M=S*fft(D*A)./sqrt(m);
end