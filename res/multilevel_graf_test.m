clear;
clc;

%%
D = 2;          % cluster diameter
X = [0 0]';     % receiver cluster center
Y = [2*D 0]';   % source cluster center

depth = 15;

l = 0;
for l = 1 : depth
    tree(l).D = D;
    tree(l).X = X;
    tree(l).Y = Y;
    
    D = D/2;
    X = X + .5 * ([D; D] .* sign(rand(2,1)-.5));
    Y = Y + .5 * ([D; D] .* sign(rand(2,1)-.5));
end

x = tree(end).X + (rand(2,1)-.5)*tree(end).D/2;     % receiver point
y = tree(end).Y + (rand(2,1)-.5)*tree(end).D/2;     % source point

figure;
hold on;
for l = 1 : length(tree)
    plot(tree(l).X(1), tree(l).X(2), 'ro', tree(l).Y(1), tree(l).Y(2), 'bo');
end
plot(x(1), x(2), 'r*', y(1), y(2), 'b*');
axis equal;


%%
kvec = logspace(1, 3, 101);
for iK = 1 : length(kvec)
    k = kvec(iK); % wave number
    lambda = 2*pi/k;
    
    for l = 1 : depth
        tree(l).Drel = tree(l).D / lambda;
        tree(l).high = tree(l).Drel > 3;
        tree(l).M = ceil(40. + 9. * tree(l).Drel);  % as in NiHu
    end
    tree(end).high = false;
    
    % p2m is always in low freq domain
    M = tree(end).M;
    nu = (-M : M).';
    rvec = tree(end).Y - y;
    [theta_y, r_y] = cart2pol(rvec(1), rvec(2));
    tree(end).multi = besselj(-nu, k*r_y) .* exp(-1i*(-nu)*theta_y);
    tree(end).Multi = padfft(tree(end).multi);
    
    % upward pass
    for l = depth-1 : -1 : 1 % to tree(l)
        M = tree(l).M;
        nu = (-M : M).';
        rvec = tree(l).Y - tree(l+1).Y;
        [Theta, R] = cart2pol(rvec(1), rvec(2));
        tree(l).m2m = besselj(-nu, k*R) .* exp(-1i*(-nu)*Theta);
        tree(l).M2M = padfft(tree(l).m2m);
        
        % compute multi in spatial domain
        multi = zeros(size(nu));
        multi(M+1 + (-tree(l+1).M : tree(l+1).M)) = tree(l+1).multi;
        multi = conv(tree(l).m2m, multi);
        multi = multi(2*M+1 + (-M:M));
        tree(l).multi = multi;
        
        % compute Multi in spectral domain
        Mfrom = tree(l+1).M;
        multi = ifft(tree(l+1).Multi);
        multi_padded = zeros(size(tree(l).M2M));
        multi_padded(1:Mfrom+1) = multi(1:Mfrom+1);
        multi_padded(end-Mfrom+1:end) = multi(end-Mfrom+1:end);
        tree(l).Multi = fft(multi_padded) .* tree(l).M2M;
    end
    
    norm(tree(1).multi - padifft(tree(1).Multi)) / norm(tree(1).multi)
    
    % m2l pass
    M = tree(1).M;
    nu = (-M : M).';
    rvec = tree(1).X - tree(1).Y;
    [Theta, R] = cart2pol(rvec(1), rvec(2));
    m2l = besselh(nu, 2, k*R) .* exp(1i*nu*(pi+Theta));
    M2L = padfft(m2l);
    
    local = conv(m2l, tree(1).multi);
    tree(1).local = local(2*M+1 + (-M:M));
    tree(1).Local = M2L .* tree(1).Multi;
    
    norm(tree(1).local - padifft(tree(1).Local)) / norm(tree(1).local)
    
    % downward pass
    for l = 1 : depth-1 % tree(l) = from cluster
        M = tree(l).M;
        nu = (-M : M).';
        rvec = tree(l+1).X - tree(l).X;
        [Theta, R] = cart2pol(rvec(1), rvec(2));
        l2l = besselj(nu, k*R) .* exp(1i*nu*(pi+Theta));
        L2L = padfft(l2l);
        
        % spatial
        local = conv(l2l, tree(l).local);
        local = local(2*M+1 + (-tree(l+1).M:tree(l+1).M));
        tree(l+1).local = local;
        
        % spectral
        Mto = tree(l+1).M;
        Loc = tree(l).Local .* L2L;
        loc = padifft(Loc);
        loc = loc(M + 1 + (-Mto : Mto));
        tree(l+1).Local = padfft(loc);
        
        norm(tree(l+1).local - padifft(tree(l+1).Local)) / norm(tree(l+1).local)
    end
    
    % l2p
    nu = (-tree(end).M : tree(end).M).';
    rvec = x - tree(end).X;
    [Theta, R] = cart2pol(rvec(1), rvec(2));
    l2p = besselj(nu, k*R) .* exp(-1i*nu*Theta);
    
    p = l2p .' * tree(end).local;
    p2 = l2p .' * padifft(tree(end).Local);

    % single level
    M = tree(1).M;
    nu = (-M : M).';
    rvec = tree(1).Y - y;
    [Theta, R] = cart2pol(rvec(1), rvec(2));
    multi = besselj(nu, k*R) .* exp(-1i*nu*Theta);

    rvec = tree(1).X - tree(1).Y;
    [Theta, R] = cart2pol(rvec(1), rvec(2));
    m2l = besselh(nu, 2, k*R) .* exp(1i*nu*(pi+Theta));
    local = conv(m2l, flipud(multi));
    local = local(2*M+1 + (-M:M));
    
    rvec = x - tree(1).X;
    [Theta, R] = cart2pol(rvec(1), rvec(2));
    l2p = besselj(nu, k*R) .* exp(-1i*nu*Theta);
    
    p1level = l2p.' * local;
    
    
    panal = besselh(0, 2, k*norm(x-y));
    err(iK) = real(log10(abs(p/panal) - 1));
    err1(iK) = real(log10(abs(p1level/panal) - 1));
    err2(iK) = real(log10(abs(p2/panal) - 1));
    Drel(iK) = tree(1).Drel;
    fprintf(1, 'Drel: %g\tM: %g\terr: %g\terr: %g\terr: %g\n', Drel(iK), tree(1).M, err(iK), err1(iK), err2(iK));
end

%%

figure;
semilogx(Drel, err, Drel, err1, Drel, err2);
xlabel('relative cluster size');
ylabel('relative errr');
legend('multilevel', 'single level', 'full spectral');
ylim([-16, 0]);
