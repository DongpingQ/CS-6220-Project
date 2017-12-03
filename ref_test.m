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

%% The underdetermined least squares system
clear
m=128;
n=16384;
l=m*4;
U=randn(m,m)+1i*randn(m,m);
U=orth(U);
V=randn(n,m)+1i*randn(n,m);
V=orth(V);
d=10.^(-6*(0:(m-1))/(m-1));
sigma=diag(d);
A=U*sigma*V';

randones=rand(1,m);
randones(randones<0.5)=-1;
randones(randones>=0.5)=1;
p=sum(bsxfun(@times,V,randones),2)/sqrt(m);
b=A*p;
epsdenom=cond(A)*norm(p);

t=zeros(1,10);
t0=zeros(1,10);
epsori=zeros(1,10);
epsrel=zeros(1,10);
for j=1:10
    tic
    x=SRFT_ls_under(A,b,l,@conj_grad,@SRFT_ls_over);
    t(j)=toc;
    epsrel(j)=norm(x-p)/epsdenom;
    tic
    [Q,R]=qr(A');
    z=R(1:m,1:m)'\b;
    p0=Q(:,1:m)*z;
    t0(j)=toc;
    epsori(j)=norm(p0-p)/epsdenom;
end
avg_ori=mean(t0);
avg_ran=mean(t);
t_rel=avg_ori/avg_ran;
max_epsrel=max(epsrel);
max_epsori=max(epsori);