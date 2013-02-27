R = 1;
nR = 5;
mesh = create_sphere_boundary(R, nR);
% mesh = quad2tria(mesh);

% First inner resonant frequency
k0 = pi;

kvec = k0 + (-0.05:.005:0.05);
error_std = nan(size(kvec));
error_bm =  nan(size(kvec));

figure;
h = semilogy(kvec, error_std, kvec, error_bm);
legend({'Error std', 'Error bm'});

for ik = 1:length(kvec);
    fprintf('%3d / %3d\n', ik, length(kvec));
    k = kvec(ik);
    if (k < 1)
        alpha = 1i;
    else
        alpha = 1i/k;
    end
    
    % Matrix assembly
    tic;
    [Hbm, Gbm] = bemHG_bm(mesh, k, alpha);
    fprintf('BM time: '); toc;
    tic;
    [H, G] = bemHG(mesh, k, 'const');
    fprintf('STD time: '); toc;
    
    % Solution
    q = ones(size(Hbm,1),1);
    
    % Analytic solution
    p_an = ones(size(q))*(-1/(1+1i*k));
    
    % Standard solution
    p_std = (H - 0.5*eye(size(H)))  \ (G*q);
    
    % Burton-Miller solution This seems to be OK
    p_bm = (Hbm - (1/2 + alpha*1i*k/2)*eye(size(Hbm))) \ ...
        (Gbm * q + (alpha/2 - 1i/(2*k)) * q);
    
    % Error calculation
    error_std(ik) = norm(p_an - p_std)/norm(p_an);
    error_bm(ik) = norm(p_an - p_bm)/norm(p_an);
    
    set(h(1), 'YData', error_std);
    set(h(2), 'YData', error_bm);
    drawnow;
end
