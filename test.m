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
sol=zeros(100,1);
acc=zeros(100,1);
for i=1:100
    A=randn(m,n);
    b=randn(m,1);
    x=SRFT_ls_under(A,b,eps,l,@srftmult_ls,@conj_grad,@randH,@SRFT_ls_over);
    [Q,R]=qr(A');
    z=R(1:m,1:m)'\b;
    p=Q(:,1:m)*z;
    sol(i)=norm(x-p);
    acc(i)=eps*norm(p);
end
plot(sol,'r.')
hold on
plot(acc,'g.')

%% Test overdetermined least squares problem
clear
m=250;
n=50;
l=4*n;
eps=1e-6;
sol=zeros(100,1);
acc=zeros(100,1);
for i=1:100
    A=randn(m,n);
    b=randn(m,1);
    x=SRFT_ls_over(A,b,eps,l,@srftmult_ls,@conj_grad,@randH);
    p=A\b;
    sol(i)=norm(A*x-b)-norm(A*p-b);
    acc(i)=eps*norm(A*p-b);
end
plot(sol,'r.')
hold on
plot(acc,'g.')
legend('Error','Bound','Location','east')

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