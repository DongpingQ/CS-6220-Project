clear

m = 32768;    n = 2.^[6,7,8,9];
l = 4*n;

for i = 1:4
    
    U = randn(m,n(i)+1) + 1i*randn(m,n(i)+1);
    U = orth(U);
    w = U(:,n(i)+1);
    U = U(:,1:n(i));
    V = randn(n(i),n(i))+1i*randn(n(i),n(i));
    V = orth(V);
    d = 10.^(-6*(0:(n(i)-1))/(n(i)-1));
    sigma = diag(d);
    A = U*sigma*V';
    
    yrand = randn(n(i),1) + 1i*randn(n(i),1);
    Ayrand = A*yrand;
    Ayrand = Ayrand/norm(Ayrand);
    Ay = sqrt(1-1e-6)*Ayrand;
    b = (1e-3).*w+Ay;
    
    for j = 1:10

        tic
        s = l(i);
        C = GaussianProjection(A, s);
        x0 = InexactLSR(A, b, s, 1);
        t0(j) = toc;
        
        tic
        [x1,P,count] = SRFT_ls_over(A,b,1.0e-4,l(i),@conj_grad);
        t1(j) = toc;

        tic
        s = l(i);
        C = CountSketch(A, s);
        x2 = InexactLSR(A, b, s, 2);
        t2(j) = toc;

        e_gau(j) = (norm(A*x0-b) - 1e-3)/(cond(A)*1e-3);
        e_srft(j) = (norm(A*x1-b) - 1e-3)/(cond(A)*1e-3);
        e_sketch(j) = (norm(A*x2-b) - 1e-3)/(cond(A)*1e-3);
    end
    
    t_gau(i) = mean(t0);
    t_srft(i) = mean(t1);
    t_sketch(i) = mean(t2);
    
    err_gau(i) = mean(e_gau);
    err_srft(i) = mean(e_srft);
    err_sketch(i) = mean(e_sketch);
end

for i = 1:4   I(i) = i;    end

figure(1)
plot(I, log10(t_gau), '-o', I, log10(t_srft), '-*', I, log10(t_sketch), '-+', 'LineWidth', 1.5)
xlabel('grid: l = 2^{[6 7 8 9]}')
ylabel('log(time)')
legend('Gaussian','SRFT','CountSketch')
title('k-SVD running time')
axis image

figure(2)
plot(I, err_gau, '-o', I,  err_srft, '-*', I, err_sketch, '-+', 'LineWidth',1.5)
xlabel('grid: l = 2^{[6 7 8 9]}')
ylabel('log(error)')
legend('Gaussian','SRFT','CountSketch')
title('max error least square solution')
axis image