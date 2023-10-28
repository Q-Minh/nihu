clear;
clc;

D = 1;          % cluster diameter
X = [0 0]';     % receiver cluster center
Y = [2*D 0]';   % source cluster center

x = X + [1;1]/2*D;     % receiver point
y = Y + [-1;-1]/2*D;     % source point

kvec = logspace(-1, 3, 101);    % wave number vector
for iK = 1 : length(kvec)
    k = kvec(iK); % wave number
    
    %
    lambda = 2*pi/k;
    Drel(iK) = D / lambda;
    M = ceil((40. + 9. * Drel(iK)));  % as in NiHu
    
    %
    nu = (-M : M).';
    
    % P2M
    rvec = Y - y;
    [theta_y, r_y] = cart2pol(rvec(1), rvec(2));
    p2m = besselj(nu, k*r_y) .* exp(-1i*nu*theta_y);
    
    % M2L
    rvec = X - Y;
    [Theta, R] = cart2pol(rvec(1), rvec(2));
    m2l = besselh(nu, 2, k*R) .* exp(1i*nu*(pi+Theta));
    
    % L2P
    rvec = x - X;
    [theta_x, r_x] = cart2pol(rvec(1), rvec(2));
    l2p = besselj(nu, k*r_x) .* exp(-1i*nu*theta_x);
    
    % Compute
    a = p2m;
    b = conv(m2l, flipud(a));
    b = b(2*M+1 + (-M:M));
    p = l2p.' * b;
    
    % Compute with FFT
    S = 2*length(p2m);
    A = fft(flipud(a), S);
    M2L = fft(m2l, S);
    B  = M2L .* A;
    bfft = ifft(B);
    bfft = bfft(2*M+1 + (-M:M));
    pfft = l2p.' * bfft;
    
    % p2l
    rvec = X - y;
    [Theta, R] = cart2pol(rvec(1), rvec(2));
    kD = k * D;
    p2l = besselh(nu, 2, k*R) .* exp(1i*nu*(pi + Theta));
    p2 = l2p.' * p2l;
    
    % Error
    panal = besselh(0, 2, k*norm(x-y));
    err(iK) = real(log10(abs(p/panal) - 1));
    errfft(iK) = real(log10(abs(pfft/panal) - 1));
    err2(iK) = real(log10(abs(p2/panal) - 1));
    
    fprintf(1, 'Drel: %g\tM: %g\terr: %g\terr: %g\terr2: %g\n', Drel(iK), M, err(iK), errfft(iK), err2(iK));
end

figure;
semilogx(Drel, err, Drel, errfft, Drel, err2);
legend('convolution', 'fft', 'p2l+l2p', 'p2m+m2p');
xlabel('relative cluster size');
ylabel('relative errr');
ylim([-16 0]);