clear all 
format compact

% This is a test.

type = 1;
m = 800;    eps = 1.0e-1;    k = 5;

n = [10 20 40 80 160];

for i = 1:length(n)
    
    if type == 1
        s = floor(m/eps^2);
        
        A = 10*rand(n(i),m);
        b = 10*rand(n(i),1);     
    else
        s = floor(m/eps^2);
        
        U = randn(m,m) + 1i*randn(m,m);
        U = orth(U);
        V = randn(n(i),m) + 1i*randn(n(i),m);
        V = orth(V');
        V = V';
        d = 10.^(-6*(0:(m-1))/(m-1)^2);
        sigma = diag(d);
        A = 10*U*sigma*V';
        
        randones=rand(1,m);
        randones(randones<0.5) = -1;
        randones(randones>=0.5) = 1;
        p = sum(bsxfun(@times,V,randones),2)/sqrt(m);
        b = A*p;
    end
    
    tol = 1.0e-2;
    x = lsqr(A, b, tol);
    
    tic
    C = GaussianProjection(A, s);
    x1 = InexactLSR(A, b, s, 1);
    t1(i) = toc;
    nF1(i) = max( abs(x - x1) );
    

    tic
    C = CountSketch(A, s);
    x2 = InexactLSR(A, b, s, 2);
    t2(i) = toc;
    nF2(i) = max( abs(x - x2) );

end


for i = 1:length(n)    I(i) = i;    end

figure(1)
plot(I, log10(t1), '-o', I, log10(t2), '-*','LineWidth',1.5)
% xlabel('grid: n = 2^{[6 7 8 9 10]}')
ylabel('log(time)')
legend('Gaussian','CountSketch')
title('least square run time')
axis image

figure(2)
plot(I, log10(nF1), I, log10(nF2), '-*','LineWidth',1.5)
% xlabel('grid: n = 2^{[6 7 8 9 10]}')
ylabel('log(error)')
legend('Gaussian','CountSketch')
title('least square errors')
axis image

