clear all 
format compact

m = 50;    eps = 1.0e-1;

s = floor(m/eps^2);

n = [200 400 800 1600 3200 6400];

for i = 1:length(n)
    
    A = 10*sin( rand(m, n(i)) );
    
    tic
    C = GaussianProjection(A, s);
    t1(i) = toc;
    
    Cinv = pinv(C, eps);

    X = A - C*Cinv*A;

    nF1(i) = norm(X,'fro');
%     fprintf('Frobenius norm: ')
%     disp(nF1(i))

    tic
    C = CountSketch(A, s);
    t2(i) = toc;

    Cinv = pinv(C, eps);

    X = A - C*Cinv*A;

    nF2(i) = norm(X,'fro');

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

