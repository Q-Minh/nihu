clear;
close all;

%% Create mesh of Helmholtz resonator

Lout = [1 1 1];
Lin = [.6 .6 .6];
Lop = .2;
Le = 1e-1;
mesh = create_brick(Lout, ceil(Lout/Le));
mesh = translate_mesh(mesh, -Lout/2);

eps = 1e-3;
mesh = mesh_section(mesh, [-1-eps; 1+eps] * Lin/2, 'nall');
mesh = mesh_section(mesh, [-Lop/2 -Lop/2 0; Lop/2 Lop/2 Inf]*(1+eps), 'nall');
mesh = drop_mesh_IDs(drop_unused_nodes(get_boundary(mesh)));
mesh = quad2tria(mesh);

[r_nodes, r_elem] = extract_core_mesh(mesh);

%% Frequency range
c = 340;
kmax = min(mesh_kmax(mesh, 8));
fmax = kmax * c / (2*pi);
fprintf('Max mesh freq: %.1f Hz\n', fmax);
fvec = 30 : 45;
nFreq = length(fvec);

%%

nElem = size(mesh.Elements, 1);

ps_conv = zeros(nElem, nFreq);
ps_bm = zeros(nElem, nFreq);
ps_hs = zeros(nElem, nFreq);

qs_ana = ones(nElem, 1);

r_c = centnorm(mesh);

for iFreq = 1 : nFreq
    
    fprintf(1, 'Freq %d from %d\n', iFreq, nFreq);
    
    f = fvec(iFreq);
    k = 2*pi*f/c;
    
    [G, H, Ht, D, ~, ~] =...
        singular_bem_3d_field(r_nodes, r_elem, [], [], k);
    
    gkernel = @(x, nx, y, ny)helmholtz_3d_slp_kernel(x, nx, y, ny, k);
    gsing = @(x, nx, corners)helmholtz_3d_slp_singular(x, nx, corners, k);
    
    hkernel = @(x, nx, y, ny)helmholtz_3d_dlp_kernel(x, nx, y, ny, k);
    hsing = @(x, nx, corners)helmholtz_3d_dlp_singular(x, nx, corners, k);
    
    htkernel = @(x, nx, y, ny)helmholtz_3d_dlpt_kernel(x, nx, y, ny, k);
    htsing = @(x, nx, corners)helmholtz_3d_dlpt_singular(x, nx, corners, k);
    
    dkernel = @(x, nx, y, ny)helmholtz_3d_hsp_kernel(x, nx, y, ny, k);
    dsing = @(x, nx, corners)helmholtz_3d_hsp_singular(x, nx, corners, k);
    
    G = bem_matrix(mesh, gkernel, gsing, G);
    H = bem_matrix(mesh, hkernel, hsing, H);
    Ht = bem_matrix(mesh, htkernel, htsing, Ht);
    D = bem_matrix(mesh, dkernel, dsing, D);
    
    %// acoustic pressure on the surface and in the field points
    ps_conv(:, iFreq) = H \ (G * qs_ana);               %// conventional
    ps_hs(:, iFreq) = D \ (Ht * qs_ana);
    coup = 1i/k;   %// coupling constant
    ps_bm(:, iFreq) = (H + coup * D) \ (G * qs_ana + coup * Ht * qs_ana);
end

%%
figure;
semilogy(fvec, err_s_conv, fvec, err_s_bm, fvec, err_s_hs, 'LineWidth', 1.2);
legend({'Conventional', 'Burton-Miller', 'Hypersingular'});
xlabel('Frequency [Hz]');
ylabel('Relative error');
%%
f_sel = 36;
f_idx = find(fvec == f_sel, 1, 'first');
figure;
subplot(2,2,1);
plot_mesh(mesh, abs(ps_conv(:,f_idx)));
axis equal;
colorbar;
title('Conventional');
subplot(2,2,2);
plot_mesh(mesh, abs(ps_bm(:,f_idx)));
axis equal;
colorbar;
title('Burton');
subplot(2,2,3);
plot_mesh(mesh, abs(ps_hs(:,f_idx)));
axis equal;
colorbar;
title('Hyper');

%%
figure;
plot(fvec, 20*log10(max(abs(ps_hs), [], 1)), ...
    fvec, 20*log10(max(abs(ps_bm), [], 1)), ...
    fvec, 20*log10(max(abs(ps_conv), [], 1)));
legend('Hypersingular', 'BM', 'conventional');