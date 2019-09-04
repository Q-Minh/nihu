clear;
clear all;

%// sphere radiator with radius R = 1, divided into 8 elements along the radius.
surface = create_slab([1, 1], [10, 10]);

%%
covariance_2d_bbfmm_matlab('init');

[r_nodes, r_elems] = extract_core_mesh(surface);
r_elems(:,1) = 22404;

%%
covariance_2d_bbfmm_matlab('mesh', r_nodes, r_elems);
%%
n_levels = 3;
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

%%
covariance_2d_bbfmm_matlab('cleanup');

%%
clear mex;