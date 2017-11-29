clear all 
format compact

m = 10;    eps = 1.0e-1;    k = 5;

s = floor(m/eps^2);

n = [400 800 1600 3200];

for i = 1:length(n)
    
    A = 10*sin( rand(m, n(i)) );
    [U, S, V] = svds(A, k);
%     A = kSVDMatrix(n(i), 4*m, 3*m, m);
    
    tic
    C = GaussianProjection(A, s);
    [U1, S1, V1] = kSVDPrototype(A, k, C);
    t1(i) = toc;
    
    nF1(i) = norm( S - S1, Inf )
    
%     Cinv = pinv(C, eps);
%     X = A - C*Cinv*A;
%     nF1(i) = norm(X,'fro');

    tic
    C = CountSketch(A, s);
    [U1, S1, V1] = kSVDPrototype(A, k, C);
    t2(i) = toc;
    
    nF2(i) = norm( S - S1, Inf )

%     Cinv = pinv(C, eps);
%     X = A - C*Cinv*A;
%     nF2(i) = norm(X,'fro');

end


for i = 1:length(n)    I(i) = i;    end

figure(1)
subplot(1,2,1)
plot(I*t1(1), t1, '-o', I*t1(1), t2, '-*')
xlabel('scaled size n')
ylabel('time')
legend('Gaussian','CountSketch')
axis image

subplot(1,2,2)
plot(I, log10(nF1), I, log10(nF2), '-*')
xlabel('size n')
ylabel('log(error)')
legend('Gaussian','CountSketch')
axis image

