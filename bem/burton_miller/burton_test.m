R = 1;
nR = 5;
mesh = quad2tria(create_sphere_boundary(R, nR));

k = 0.75;
alpha = 0;

%%
[Hbm, Gbm] = bemHG_bm(mesh, k, alpha);
[H, G] = bemHG(mesh, k, 'const');

%%
q = ones(size(Hbm,1),1);
p_bm = (Hbm - (1/2 + 1i*k/2*alpha)*eye(size(Hbm))) \ ...
       (Gbm * q + (alpha/2 - 1i/2/k) * q);

p = (H - 0.5*eye(size(H)))  \ (G*q);   
pbm = (Hbm - 0.5*eye(size(Hbm))) \ (Gbm * q - 1i/(2*k)*q);
%%
subplot(1,2,1);
plot_mesh(mesh, real(pbm)); colorbar;
subplot(1,2,2);
plot_mesh(mesh, real(p)); colorbar;
