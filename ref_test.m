% This test file tries to imitate the examples in the referenced papers
%% The overdetermined least squares system
clear
m=65536;
n=256;
l=n*4;
U=randn(m,n+1)+1i*randn(m,n+1);
U=orth(U);
w=U(:,n+1);
U=U(:,1:n);
V=randn(n,n)+1i*randn(n,n);
V=orth(V);
d=10.^(-6*(0:(n-1))/(n-1));
sigma=diag(d);
A=U*sigma*V';
yrand=randn(n,1)+1i*randn(n,1);
Ayrand=A*yrand;
Ayrand=Ayrand/norm(Ayrand);
Ay=sqrt(1-1e-6)*Ayrand;
b=(1e-3).*w+Ay;

t=zeros(1,10);
t0=zeros(1,10);
condnum=zeros(1,10);
iternum=zeros(1,10);
epsrel=zeros(1,10);
for j=1:10
    tic
    [x,P,count]=SRFT_ls_over(A,b,0.5e-10,l,@conj_grad);
    t(j)=toc;
    tic
    y=A\b;
    t0(j)=toc;
    condnum(j)=cond(A/P);
    iternum(j)=count;
    del=norm(A*x-b);
    epsrel(j)=(del-1e-3)/(cond(A)*1e-3);
end
avg_ori=mean(t0);
avg_ran=mean(t);
t_rel=avg_ori/avg_ran;
max_cond=max(condnum);
max_iter=max(iternum);
max_rel=max(epsrel);
