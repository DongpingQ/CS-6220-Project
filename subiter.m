function Q = subiter(A,q,l)
[~,n]=size(A);
omega=normrnd(0,1,n,l);
Y0=A*omega;
[Q0,~]=qr(Y0);
for j=1:q
    Yhat=A'*Q0;
    [Qhat,~]=qr(Yhat);
    Y=A*Qhat;
    [Q1,~]=qr(Y);
    Q0=Q1;
end
Q=Q1;
end