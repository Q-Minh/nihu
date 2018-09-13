clear;
L = [1 .1 .1];
Le = .1/2;
mesh = quad2tria(create_brick_boundary(L, L/Le));
freq = 100;
c = 340;
om = 2*pi*freq;
k = om / c;

gkernel = @(x, nx, y, ny)helmholtz_3d_slp_kernel(x, nx, y, ny, k);
gsing = @(x, nx, corners)helmholtz_3d_slp_singular(x, nx, corners, k);

hkernel = @(x, nx, y, ny)helmholtz_3d_dlp_kernel(x, nx, y, ny, k);
hsing = @(x, nx, corners)helmholtz_3d_dlp_singular(x, nx, corners, k);


G = bem_matrix(mesh, gkernel, gsing);
H = bem_matrix(mesh, hkernel, hsing);
