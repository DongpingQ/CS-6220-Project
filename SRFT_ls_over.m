function x = SRFT_ls_over(A,b,eps,l,conj_grad)
% This algorithm solves overdetermined Least Squares Problem Ax=b with
% A of size m by n and m >= n by means of SRFT
% Input: A: target matrix
%        b: target result vector
%        eps: tolerance for conjugate gradient iterations
%        l: parameter used in the algorithm m >= l > n
%        srftmultback, conj_grad, randH: helper functions
% Output: x: the solution to the overdetermined system

%%
% Set-up and construct random matrix
[m,~]=size(A);

S=randi(m,1,l);
D=2*pi*rand(m,1);
D=cos(D)+1i*sin(D);
% Z1=2*pi*rand(m,1);
% Z2=2*pi*rand(m,1);
% Z1=cos(Z1)+1i*sin(Z1);
% Z2=cos(Z2)+1i*sin(Z2);
% Pi1=randperm(m);
% Pi2=randperm(m);
% 
% Y=bsxfun(@times,A,Z2);
% z=Z2.*b;
% Y=Y(Pi2,:);
% z=z(Pi2,:);
% tic
% for j=1:m-1
%     th=2*pi*rand;
%     givens=[cos(th) sin(th);-sin(th) cos(th)];
%     Y(m-j:m-j+1,:)=givens*Y(m-j:m-j+1,:);
%     z(m-j:m-j+1,:)=givens*z(m-j:m-j+1,:);
% end
% toc
% Y=bsxfun(@times,Y,Z1);
% z=Z1.*z;
% Y=Y(Pi1,:);
% z=z(Pi1,:);
% for j=1:m-1
%     th=2*pi*rand;
%     givens=[cos(th) sin(th);-sin(th) cos(th)];
%     Y((m-j):(m-j+1),:)=givens*Y((m-j):(m-j+1),:);
%     z((m-j):(m-j+1),:)=givens*z((m-j):(m-j+1),:);
% end
% E=fft(bsxfun(@times,Y,D));
% z=fft(D.*z);
E=fft(bsxfun(@times,A,D));
E=E(S,:);
E=E/sqrt(m);
z=fft(D.*b);
z=z(S,:);
z=z/sqrt(m);


%%
% Least squares problem algorithm
[Q,R,Per]=qr(E,0);
P=R(:,Per);
z=P\(Q'*z);
x=conj_grad(A,P,b,z,eps);
end