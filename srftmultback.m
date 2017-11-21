function [R,M] = srftmultback(l,A)
% This function generates a l*m sized subsampled random Fourier transform
% R=SFD and computes the matrix multiplication RA

[m,~]=size(A);
% Generate S
S = zeros(l,m);
for j = 1:l
    sj=randi(m);
    S(j,sj)=1;
end

% Generate D
d=zeros(1,m);
for j = 1:m
    theta=2*pi*rand;
    d(j)=cos(theta)+sin(theta)*1i;
end
D=diag(d);

R=fft(S)*D;
M=S*fft(D*A);
end