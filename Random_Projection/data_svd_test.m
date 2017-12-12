load('A1'); load('A2'); load('A3');
load('A4'); load('A5'); load('A6');

eps = 1.0e-1;

k = 10;

for i = 1:6
    
    switch i
        case 1 
            A = A1;
        case 2 
            A = A2;
        case 3 
            A = A3;    
        case 4 
            A = A4;
        case 5 
            A = A5;
        case 6 
            A = A6;
    end
    m = length(A);
    s = floor(m/eps^2);
    
    [U, S, V] = svds(A, k);

    tic
    C = GaussianProjection(A, s);
    [U1, S1, V1] = kSVDPrototype(A, k, C);
    t1(i) = toc;

    nF1(i) = norm( S - S1, Inf );

    Cinv = pinv(C, eps);
    X = A - C*Cinv*A;
    nF11(i) = norm(X,'fro');

    tic
    C = CountSketch(A, s);
    [U1, S1, V1] = kSVDPrototype(A, k, C);
    t2(i) = toc;

    nF2(i) = norm( S - S1, Inf );

    Cinv = pinv(C, eps);
    X = A - C*Cinv*A;
    nF22(i) = norm(X,'fro');
    
    e_gau(i) = norm(u0*s0*v0'-A);
    e_srft(i) = norm(u1*s1*v1'-A);
    e_sketch(i) = norm(u2*s2*v2'-A);
end

for i = 1:6    I(i) = i;    end

figure(1)
plot(I, log10(t1), '-o', I, log10(t2), '-*','LineWidth', 1.5)
% xlabel('n = [0.4 0.8 1.6 3.2 6.4]*100')
ylabel('time')
legend('Gaussian','CountSketch')
title('k-SVD running time')
axis image

figure(2)
plot(I, log10(nF1), '-o', I,  log10(nF2), '-*','LineWidth',1.5)
% xlabel('n = [0.4 0.8 1.6 3.2 6.4]*100')
ylabel('log(error)')
legend('Gaussian','CountSketch')
title('k-SVD log(error)')
axis image

figure(3)
plot(I, log10(nF11), '-o', I, log10(nF22), '-*','LineWidth',1.5)
% xlabel('n = [0.4 0.8 1.6 3.2 6.4]*100')
ylabel('log(error)')
legend('Gaussian','CountSketch')
title('log(Frobenius norm)')
axis image