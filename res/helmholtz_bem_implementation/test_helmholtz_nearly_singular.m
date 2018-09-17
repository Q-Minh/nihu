clear;
close all;

%% Create mesh of Helmholtz resonator
mesh = quad2tria(create_sphere_boundary(1, 6));
k = 1;

K = @(x, nx, y, ny)helmholtz_3d_dlpt_kernel(x, nx, y, ny, k);
Kns = @(x, nx, corners)helmholtz_3d_dlpt_nearly_singular(x, nx, corners, k);

K = bem_matrix(mesh, K, [], Kns);
