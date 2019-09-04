clear;
clear all;

%// sphere radiator with radius R = 1, divided into 8 elements along the radius.
surface = create_slab([1, 1], [100, 100]);

%%
covariance_2d_bbfmm_matlab('init');

[r_nodes, r_elems] = extract_core_mesh(surface);
r_elems(:,1) = 22404;

%%
covariance_2d_bbfmm_matlab('mesh', r_nodes, r_elems);
%%
n_levels = 6;
covariance_2d_bbfmm_matlab('tree', n_levels);

%%
tic;
covariance_2d_bbfmm_matlab('matrix');
toc;

%%
xct = ones(size(surface.Elements, 1), 1);
tic;
res = covariance_2d_bbfmm_matlab('mvp', xct);
toc;

%% Eigenvalues

[~, ~, w] = geo2gauss(surface, 1);
W = spdiags(w, 0, size(w,1), size(w,1));

fprintf('Eigenvalues ...\n');
tic;
n_modes = 10;
Afun = @(x)covariance_2d_bbfmm_matlab('mvp', x);
[phi, lam] = eigs(Afun, size(w,1), W, n_modes, 'lm');
toc;
%%
covariance_2d_bbfmm_matlab('cleanup');

%%
%clear mex;