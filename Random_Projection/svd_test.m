clear all 
format compact

type = 2;

m = 100;    eps = 1.0e-1;    k = 5;

n = 5*[64 128 256 512 1024];

    
for i = 1:length(n)
    
    if type == 1
        s = floor(8*i);
        
        A = 10*rand(m, n(i));
    else
        s = floor(8*i);
        
        m=1024;
        n=1024;
        k=504;
        l=k+8;
        sigma = 10.^(-12*(0:l+1)'/(l+1));
        U = orth(randn(n(i),l+2) + 1i*randn(n(i),l+2));
        V = orth(randn(n(i),l+2) + 1i*randn(n(i),l+2));
        A = 10*U*diag(sigma)*V';
    end
    
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
plot(I, log10(t1), '-o', I, log10(t2), '-*','LineWidth', 1.5)
xlabel('grid: n = 2^{[6 7 8 9 10]}')
ylabel('log(time)')
legend('Gaussian','CountSketch')
title('k-SVD running time')
axis image

figure(2)
plot(I, log10(nF1), '-o', I,  log10(nF2), '-*','LineWidth',1.5)
xlabel('grid: n = 2^{[6 7 8 9 10]}')
ylabel('log(error)')
legend('Gaussian','CountSketch')
title('max error of top k singular values')
axis image

figure(3)
plot(I, log10(nF11), '-o', I, log10(nF22), '-*','LineWidth',1.5)
xlabel('grid: n = 2^{[6 7 8 9 10]}')
ylabel('log(error)')
legend('Gaussian','CountSketch')
title('Frobenius norm error')
axis image