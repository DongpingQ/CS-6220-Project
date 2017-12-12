% This test file tries to imitate the examples in the referenced papers
%% The overdetermined least squares system
clear
m=32768;
n=64;
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
    [x,P,count]=SRFT_ls_over(A,b,0.5e-14,l,@conj_grad);
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
m=256;
n=32768;
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

%% The 3rd SRFT SVD test on reference
clear
m=1024;
n=1024;
k=504;
l=k+8;
sigma=10.^(-12*(0:l+1)'/(l+1));
U=orth(randn(m,l+2)+1i*randn(m,l+2));
V=orth(randn(n,l+2)+1i*randn(n,l+2));
A=U*diag(sigma)*V';

t0=zeros(1,10);
t=zeros(1,10);
deldir=zeros(1,10);
delfas=zeros(1,10);

for j=1:10
    tic
    [u1,s1,v1]=svds(A,k);
    t0(j)=toc;
    tic
    [u2,s2,v2]=SRFT_svd(A,k);
    t(j)=toc;
    deldir(j)=norm(u1*s1*v1'-A);
    delfas(j)=norm(u2*s2*v2'-A);
end
avg_ori=mean(t0);
avg_alg=mean(t);
t_rel=avg_ori/avg_alg;
max_dir=max(deldir);
max_fas=max(delfas);

%% First One-pass SVD test from reference
clear
m=1000;
n=1000;
r=5;
gamma=0.01;
A=[eye(r) zeros(r,n-r);zeros(m-r,r) zeros(m-r,n-r)];
G=zeros(m,n);
for j=1:m
    for k=1:n
        ran=mvnrnd([0,0],eye(2));
        G(j,k)=ran(1)+1i*ran(2);
    end
end
A=A+sqrt(gamma*r/m/n)*G;
[u0,s0,v0]=svds(A,r);
relden=norm(A-u0*s0*v0','fro');

kupper=floor((r+log(n))*log(r));
lupper=floor((kupper+log(m))*log(kupper));
relerr1=zeros(kupper-r,lupper-r);
relerr2=zeros(kupper-r,lupper-r);
t1=zeros(kupper-r,lupper-r);
t2=zeros(kupper-r,lupper-r);
for k=r+1:kupper
    for l=r+1:lupper
        tic
        [u1,s1,v1]=fixedrankapprox(A,r,k,l,'G',@lowrankapprox,@Gaussian_sketch,@SRFT_sketch);
        t1(k-r,l-r)=toc;
        relerr1(k-r,l-r)=norm(A-u1*s1*v1','fro')/relden-1;
        tic
        [u2,s2,v2]=fixedrankapprox(A,r,k,l,'S',@lowrankapprox,@Gaussian_sketch,@SRFT_sketch);
        t2(k-r,l-r)=toc;
        relerr2(k-r,l-r)=norm(A-u2*s2*v2','fro')/relden-1;
    end
end
mesh(t1-t2)
figure
mesh(relerr1,'edgecolor','b')
hold on
mesh(relerr2,'edgecolor','r')
alpha(0.3)

%% Second one-pass test case from reference
n=1000;
p=2;
A=diag((1:n).^(-p));
r=5;
[u0,s0,v0]=svds(A,r);
relden=norm(A-u0*s0*v0','fro');

kupper=floor((r+log(n))*log(r));
lupper=floor((kupper+log(n))*log(kupper));
relerr1=zeros(kupper-r,lupper-r);
relerr2=zeros(kupper-r,lupper-r);
t1=zeros(kupper-r,lupper-r);
t2=zeros(kupper-r,lupper-r);
for k=r+1:kupper
    for l=r+1:lupper
        tic
        [u1,s1,v1]=fixedrankapprox(A,r,k,l,'G',@lowrankapprox,@Gaussian_sketch,@SRFT_sketch);
        t1(k-r,l-r)=toc;
        relerr1(k-r,l-r)=norm(A-u1*s1*v1','fro')/relden-1;
        tic
        [u2,s2,v2]=fixedrankapprox(A,r,k,l,'S',@lowrankapprox,@Gaussian_sketch,@SRFT_sketch);
        t2(k-r,l-r)=toc;
        relerr2(k-r,l-r)=norm(A-u2*s2*v2','fro')/relden-1;
    end
end
mesh(t1-t2)
figure
mesh(relerr1,'edgecolor','b')
hold on
mesh(relerr2,'edgecolor','r')
alpha(0.3)