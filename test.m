%% Test SVD through SRFT
m=20;
n=30;
k=5;
err=zeros(100,1);
sv=zeros(100,1);
for i=1:100
    A=randn(m,n);
    [U,S,V]=SRFT_svd(A,k,@srftmultback);
    err(i)=norm(U*S*V'-A);
    s=svds(A,k+1);
    %sv(i)=((sqrt(k*(n-k)+1)+1)*sqrt(2*m+1)+sqrt(k*(n-k)+1)*sqrt(2*m))*s(k+1);
    %sv(i)=sqrt(k*m*n)*s(k+1);
    sv(i)=sqrt(max(m,n))*s(k+1);
end
plot(err,'r.')
hold on
plot(sv,'g.')

%% Test underdetermined least squares problem
clear
m=5;
n=25;
l=4*m;
eps=1e-6;
sol1=zeros(100,1);
acc=zeros(100,1);
for i=1:100
    A=randn(m,n);
    b=randn(m,1);
    x=SRFT_ls_under(A,b,eps,l,@srftmult_ls,@conj_grad,@randH,@SRFT_ls_over);
    [Q,R]=qr(A');
    z=R(1:m,1:m)'\b;
    p=Q(:,1:m)*z;
    sol1(i)=norm(x-p);
    acc(i)=eps*norm(p);
end
plot(sol1,'r.')
hold on
plot(acc,'g.')

%% Test overdetermined least squares problem
clear
m=250;
n=50;
l=4*n;
eps=1e-6;
sol1=zeros(100,1);
acc=zeros(100,1);
sol2=zeros(100,1);
t1=zeros(100,1);
t2=zeros(100,1);
for i=1:100
    A=randn(m,n);
    b=randn(m,1);
    tic
    x=SRFT_ls_over(A,b,eps,l,@conj_grad,@randH);
    t1(i)=toc;
    tic
    y=SRFT_ls_overold(A,b,eps,l,@srftmult_ls,@conj_grad,@randH);
    t2(i)=toc;
    p=A\b;
    sol1(i)=norm(A*x-b)-norm(A*p-b);
    acc(i)=eps*norm(A*p-b);
    sol2(i)=norm(A*y-b)-norm(A*p-b);
end
plot(sol1,'r.')
hold on
plot(sol2,'b.')
hold on
plot(acc,'g.')
legend('Error1','Error2','Bound','Location','east')
mean(t1)
mean(t2)

%% Another large-scale overdetermined system
clear
m=1000;
n=100;
l=n*4;
eps=1e-6;
A=randn(m,n);
b=randn(m,1);
tic
x=SRFT_ls_over(A,b,eps,l,@conj_grad,@randH);
toc
tic
y=SRFT_ls_overold(A,b,eps,l,@srftmult_ls,@conj_grad,@randH);
toc
% p=A\b;
% norm(A*x-b)-norm(A*p-b)
% norm(A*y-b)-norm(A*p-b)

%% Test one-pass SVD with Gaussian and SRFT
m=150;
n=200;
r=30;
k=2*r+1;
l=4*r+2;
err1=zeros(100,1);
err2=zeros(100,1);
elp1=zeros(100,1);
elp2=zeros(100,1);
for i=1:100
    A=randn(m,n);
    tic;
    [Q1,S1,V1]=fixedrankapprox(A,r,k,l,'G',@lowrankapprox,@Gaussian_sketch,@SRFT_sketch,@srftmultfront,@srftmultback);
    elp1(i)=toc;
    tic
    [Q2,S2,V2]=fixedrankapprox(A,r,k,l,'S',@lowrankapprox,@Gaussian_sketch,@SRFT_sketch,@srftmultfront,@srftmultback);
    elp2(i)=toc;
    err1(i)=norm(Q1*S1*V1'-A,'fro');
    err2(i)=norm(Q2*S2*V2'-A,'fro');
end
mean(elp1)
mean(elp2)
plot(err1,'r.')
hold on
plot(err2,'g.')