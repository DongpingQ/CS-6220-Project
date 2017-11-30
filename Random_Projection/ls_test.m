clear all 
format compact

m = 10;    eps = 1.0e-1;    k = 5;    tol = eps;

s = floor(m/eps^2);

n = [400 800 1600 3200 6400];

for i = 1:length(n)
    
    A = 10*rand(n(i),m);
    b = 10*rand(n(i),1);  
    x = lsqr(A, b, tol);
%     A = kSVDMatrix(n(i), 4*m, 3*m, m);
    
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
plot(I*t1(1), t1, '-o', I*t1(1), t2, '-*','LineWidth',1.5)
xlabel('scaled size n')
ylabel('time')
legend('Gaussian','CountSketch')
axis image

figure(2)
plot(I, log10(nF1), I, log10(nF2), '-*','LineWidth',1.5)
xlabel('size n')
ylabel('log(error)')
legend('Gaussian','CountSketch')
axis image

