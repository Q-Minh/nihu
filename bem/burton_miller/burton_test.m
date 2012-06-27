R = 1;
nR = 6;
mesh = quad2tria(create_sphere_boundary(R, nR));

k = 0.5;
alpha = 1i/k;

%%
[Hbm, Gbm] = bemHG_bm(mesh, k, alpha);
[H, G] = bemHG(mesh, k, 'const');

%%
q = ones(size(Hbm,1),1);
p_bm = (Hbm - (1 + 1i*k/2*alpha)*eye(size(Hbm))) \ ...
       (Gbm * q + (alpha - 1i/2/k) * q);

p = (H - eye(size(H)))  \ (G*q);   
%%
subplot(1,2,1);
plot_mesh(mesh, abs(p_bm)); colorbar;
subplot(1,2,2);
plot_mesh(mesh, abs(p)); colorbar;
