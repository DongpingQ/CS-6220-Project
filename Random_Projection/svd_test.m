clear all 
format compact

m = 10;    eps = 1.0e-1;    k = 5;

s = floor(m/eps^2);

n = [400 800 1600 3200 6400];

for i = 1:length(n)
    
    A = 10*rand(m, n(i));
    [U, S, V] = svds(A, i*k);
    
%     A = kSVDMatrix(n(i), 4*m, 3*m, m);
    
    tic
    C = GaussianProjection(A, s);
    [U1, S1, V1] = kSVDPrototype(A, i*k, C);
    t1(i) = toc;
    
    nF1(i) = norm( S - S1, Inf );
    
    Cinv = pinv(C, eps);
    X = A - C*Cinv*A;
    nF11(i) = norm(X,'fro');

    tic
    C = CountSketch(A, s);
    [U1, S1, V1] = kSVDPrototype(A, i*k, C);
    t2(i) = toc;
    
    nF2(i) = norm( S - S1, Inf );

    Cinv = pinv(C, eps);
    X = A - C*Cinv*A;
    nF22(i) = norm(X,'fro');

end


for i = 1:length(n)    I(i) = i;    end

figure(1)
plot(I*t1(1), t1, '-o', I*t1(1), t2, '-*','LineWidth', 1.5)
xlabel('n = [0.4 0.8 1.6 3.2 6.4]*100')
ylabel('time')
legend('Gaussian','CountSketch')
title('k-SVD running time')
axis image

figure(2)
plot(I, log10(nF1), '-o', I,  log10(nF2), '-*','LineWidth',1.5)
xlabel('n = [0.4 0.8 1.6 3.2 6.4]*100')
ylabel('log(error)')
legend('Gaussian','CountSketch')
title('k-SVD log(error)')
axis image

figure(3)
plot(I, log10(nF11), '-o', I, log10(nF22), '-*','LineWidth',1.5)
xlabel('n = [0.4 0.8 1.6 3.2 6.4]*100')
ylabel('log(error)')
legend('Gaussian','CountSketch')
title('log(Frobenius norm)')
axis image