clear 

m = 1024;    n = 1024;
k = 2.^[4,5,6,7,8,9] - 8;
l = k + 8;

eps = 1.0e-1;

for i = 1:6
    
    sigma = 10.^(-12*(0:l(i)+1)'/(l(i)+1));
    U = orth(randn(m,l(i)+2) + 1i*randn(m,l(i)+2));
    V = orth(randn(n,l(i)+2) + 1i*randn(n,l(i)+2));
    A = U*diag(sigma)*V';
    
    for j = 1:10
        tic
        s = l(i);
        Q = GaussianProjection(A, s);
        [u0,s0,v0] = kSVDPrototype(A, k(i), Q);
        t0(j) = toc;
        
        tic
        [u1,s1,v1] = SRFT_svd(A,k(i));
        t1(j) = toc;
        
        tic
        s = k(i);
        C = CountSketch(A, s);
        [u2,s2,v2] = kSVDPrototype(A, k(i), C);
        t2(j) = toc;
        
        e_gau(j) = norm(u0*s0*v0'-A);
        e_srft(j) = norm(u1*s1*v1'-A);
        e_sketch(j) = norm(u2*s2*v2'-A);
    end

    t_gau(i) = mean(t0);
    t_srft(i) = mean(t1);
    t_sketch(i) = mean(t2);
    err_gau(i) = max(e_gau);
    err_srft(i) = max(e_srft);
    err_sketch(i) = max(e_sketch);
end

for i = 1:6    I(i) = i;    end

figure(1)
plot(I, log10(t_gau), '-o', I, log10(t_srft), '-*', I, log10(t_sketch), '-+', 'LineWidth', 1.5)
xlabel('grid: l = 2^{[4 5 6 7 8 9]}')
ylabel('log(time)')
legend('Gaussian','SRFT','CountSketch')
title('k-SVD running time')
axis image

figure(2)
plot(I, log10(err_gau), '-o', I,  log10(err_srft), '-*', I, log10(err_sketch), '-+', 'LineWidth',1.5)
xlabel('grid: l = 2^{[4 5 6 7 8 9]}')
ylabel('log(error)')
legend('Gaussian','SRFT','CountSketch')
title('max error of singular values')
axis image