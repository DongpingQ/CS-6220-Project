function [Y,W,psi] = SRFT_sketch(A,k,l)
% This function performs a sketch of A with SRFT

Y=srftmultfront(k,A);
[psi,W]=srftmultback(l,A);
end