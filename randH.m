function H = randH(m)
% This function constructs another random matrix H to use together with SRFT
% H = th1*perm1*z1*th2*perm2*z2 and compute T = SRFT*H
perm1=eye(m);
perm2=eye(m);
idx1=randperm(m);
idx2=randperm(m);
perm1=perm1(idx1,:);
perm2=perm2(idx2,:);

c1=zeros(1,m);
c2=zeros(1,m);
for j = 1:m
    theta1=2*pi*rand;
    theta2=2*pi*rand;
    c1(j)=cos(theta1)+sin(theta1)*1i;
    c2(j)=cos(theta2)+sin(theta2)*1i;
end
z1=diag(c1);
z2=diag(c2);

th1=eye(m);
th2=eye(m);
for j = 1:m-1
    t1=2*pi*rand;
    t2=2*pi*rand;
    temp1=eye(m);
    temp2=eye(m);
    temp1(j,j)=cos(t1);
    temp1(j,j+1)=sin(t1);
    temp1(j+1,j)=-sin(t1);
    temp1(j+1,j+1)=cos(t1);
    temp2(j,j)=cos(t2);
    temp2(j,j+1)=sin(t2);
    temp2(j+1,j)=-sin(t2);
    temp2(j+1,j+1)=cos(t2);
    th1=th1*temp1;
    th2=th2*temp2;
end
H=th1*perm1*z1*th2*perm2*z2;
end