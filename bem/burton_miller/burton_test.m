R = 1;
nR = 3;
mesh = create_sphere_boundary(R, nR);

k = 1;

[H, G] = bemHG_bm(mesh, k);

figure;
subplot(1,2,1);  pcolor(abs(H)); shading flat; axis equal;
subplot(1,2,2);  pcolor(abs(G)); shading flat; axis equal;