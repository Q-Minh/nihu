clear;
close all;

%% Create mesh of Helmholtz resonator
mesh = quad2tria(create_slab([1, 1], [10 10]));

D = @(x, nx, y, ny)laplace_3d_hsp_kernel(x, nx, y, ny);
Dns = @(x, nx, corners)laplace_3d_hsp_nearly_singular(x, nx, corners);

D = bem_matrix(mesh, D, [], Dns);
