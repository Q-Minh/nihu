clear;
close all;

Lout = [.2 .2 1];
Le = 5e-2/1;
eps = 1e-3;

tria = true;

mesh = create_brick(Lout, ceil(Lout/Le));
mesh = translate_mesh(mesh, -Lout/2);

mesh = drop_mesh_IDs(drop_unused_nodes(get_boundary(mesh)));
mesh = flip_mesh(mesh);

if (tria)
    mesh = quad2tria(mesh);
end

c = 340;
kmax = min(mesh_kmax(mesh, 8));

fmax = kmax * c / (2*pi);
fprintf('Max mesh freq: %.1f Hz\n', fmax);

[~, exc_sel] = mesh_select(mesh, sprintf('abs(z) > %f', Lout(3)/2-eps), 'ind', 'all');
%%
%// call C++ code at wave number k \approx pi
fvec = 100;
nFreq = length(fvec);

nElem = size(mesh.Elements, 1);

ps_conv = zeros(nElem, nFreq);

ps_bm = zeros(nElem, nFreq);

ps_hs = zeros(nElem, nFreq);

qs_ana = zeros(nElem, 1);
qs_ana(exc_sel) = 1;

ps_ana = zeros(nElem, nFreq);

r_c = centnorm(mesh);

err_s_conv = zeros(1, nFreq);

err_s_hs = zeros(1, nFreq);
err_s_bm = zeros(1, nFreq);


for iFreq = 1 : nFreq
    
    progbar(1, nFreq, iFreq);
    f = fvec(iFreq);
    k = 2*pi*f/c;
    
    gkernel = @(x, nx, y, ny)helmholtz_3d_slp_kernel(x, nx, y, ny, k);
    gsing = @(x, nx, corners)helmholtz_3d_slp_singular(x, nx, corners, k);

    hkernel = @(x, nx, y, ny)helmholtz_3d_dlp_kernel(x, nx, y, ny, k);
    hsing = @(x, nx, corners)helmholtz_3d_dlp_singular(x, nx, corners, k);

    htkernel = @(x, nx, y, ny)helmholtz_3d_dlpt_kernel(x, nx, y, ny, k);
    htsing = @(x, nx, corners)helmholtz_3d_dlpt_singular(x, nx, corners, k);
    
    dkernel = @(x, nx, y, ny)helmholtz_3d_hsp_kernel(x, nx, y, ny, k);
    dsing = @(x, nx, corners)helmholtz_3d_hsp_singular(x, nx, corners, k);

    G = bem_matrix(mesh, gkernel, gsing);
    H = bem_matrix(mesh, hkernel, hsing);
    Ht = bem_matrix(mesh, htkernel, htsing);
    D = bem_matrix(mesh, dkernel, dsing);
    
    H = H - 0.5*eye(nElem);
    Ht = Ht + 0.5*eye(nElem);
    
    %// acoustic pressure on the surface and in the field points
    ps_conv(:, iFreq) = H \ (G * qs_ana);               %// conventional
    coup = 1i/k;   %// coupling constant
    ps_bm(:, iFreq) = (H + coup * D) \ (G * qs_ana + coup * Ht * qs_ana);
    
    ps_hs(:, iFreq) = (D) \ (Ht * qs_ana);
    
    % Solve analytically
    A = [ 1i*k*exp(-1i*k*Lout(3)/2), -(1i*k)*exp(1i*k*Lout(3)/2);
        -1i*k*exp(-1i*k*(-Lout(3)/2)), (1i*k)*exp(1i*k*(-Lout(3)/2)) ];
    sol = A \ [1; 1];
    pp = sol(1); pm = sol(2);
    
    ps_ana(:, iFreq) = pp*exp(-1i*k*r_c(:,3)) + pm*exp(1i*k*r_c(:,3));
    
    err_s_conv(:, iFreq) = norm(ps_conv(:, iFreq) - ps_ana(:, iFreq)) / norm(ps_ana(:, iFreq));
    err_s_hs(:, iFreq) = norm(ps_hs(:, iFreq) - ps_ana(:, iFreq)) / norm(ps_ana(:, iFreq));
    err_s_bm(:, iFreq) = norm(ps_bm(:, iFreq) - ps_ana(:, iFreq)) / norm(ps_ana(:, iFreq));
    
    
end
%%
if (0)
    figure;
    subplot(3, 1, 1);
    plot(fvec, real(pf_conv(1, :)), ...
        fvec, real(pf_bm(1, :)), ...
        fvec, real(pf_hs(1, :)));
    hold on;
    plot(fvec, real(pf_ana(1,:)), 'k');
    subplot(3, 1, 2);
    plot(fvec, imag(pf_conv(1, :)), ...
        fvec, imag(pf_bm(1, :)), ...
        fvec, imag(pf_hs(1, :)));
    hold on;
    plot(fvec, imag(pf_ana(1,:)), 'k');
    subplot(3, 1, 3);
    plot(fvec, log10(abs(pf_conv(1, :))), ...
        fvec, log10(abs(pf_bm(1, :))), ...
        fvec, log10(abs(pf_hs(1, :))));
    hold on;
    plot(fvec, log10(abs(pf_ana(1,:))), 'k');
    % // plot results
end
%%
figure;
semilogy(fvec, err_s_conv, fvec, err_s_bm, fvec, err_s_hs, 'LineWidth', 1.2);
legend({'Conventional', 'Burton-Miller', 'Hypersingular'});
xlabel('Frequency [Hz]');
ylabel('Relative error');
%%
if (1)
    f_sel = 100;
    f_idx = find(fvec == f_sel, 1, 'first');
    figure;
    subplot(2,2,1);
    plot_mesh(mesh, real(ps_conv(:,f_idx)));
    axis equal;
    colorbar;
    caxis([-1 1]*max(abs(ps_ana(:, f_idx))));
    title('Conventional');
    subplot(2,2,2);
    plot_mesh(mesh, real(ps_bm(:,f_idx)));
    axis equal;
    colorbar;
    caxis([-1 1]*max(abs(ps_ana(:, f_idx))));
    title('Burton');
    subplot(2,2,3);
    plot_mesh(mesh, real(ps_hs(:,f_idx)));
    axis equal;
    colorbar;
    caxis([-1 1]*max(abs(ps_ana(:, f_idx))));
    title('Hyper');
    subplot(2,2,4);
    plot_mesh(mesh, real(ps_ana(:,f_idx)));
    axis equal;
    colorbar;
    caxis([-1 1]*max(abs(ps_ana(:, f_idx))));
    title(sprintf('Analytical f = %g Hz', f_sel));
end

%%
figure
plot(r_c(:, 3), real(ps_conv(:, f_idx)), '*'); hold on;
plot(r_c(:, 3), real(ps_ana(:, f_idx)), 'r*');
grid on;