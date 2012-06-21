R = 1;
nR = 3;
mesh = create_sphere_boundary(R, nR);

k = 1;

[H, G] = bemHG_bm(mesh, k);