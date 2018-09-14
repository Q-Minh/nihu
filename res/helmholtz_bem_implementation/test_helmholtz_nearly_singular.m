clear;
close all;

%% Create mesh of Helmholtz resonator
mesh = quad2tria(create_slab([1, 1], [10 10]));
k = 1;

D = @(x, nx, y, ny)helmholtz_3d_hsp_kernel(x, nx, y, ny, k);
Dns = @(x, nx, corners)helmholtz_3d_hsp_nearly_singular(x, nx, corners, k);

D = bem_matrix(mesh, D, [], Dns);
