R = 1;
nR = 7;
mesh = create_sphere_boundary(R, nR);

k = 1;
alpha = 1i/k;

%%
[Hbm, Gbm] = bemHG_bm(mesh, k, alpha);
[H, G] = bemHG(mesh, k, 'lin');

%%
q = ones(size(H,1),1);
p = (H - .5*eye(size(H))) \ G * q;
p_bm = (Hbm - .5*eye(size(Hbm))) \ (Gbm * q + alpha/2 * q);

%%
subplot(1,2,1);
plot_mesh(mesh, abs(p_bm)); colorbar;
subplot(1,2,2);
plot_mesh(mesh, abs(p)); colorbar;
