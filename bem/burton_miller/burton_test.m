R = 1;
nR = 5;
mesh = quad2tria(create_sphere_boundary(R, nR));

k = 0.7*min(bemkmax(mesh, 7));

%k = 1;
alpha = 1i/k ;

%%
[Hbm, Gbm] = bemHG_bm(mesh, k, alpha);
[H, G] = bemHG(mesh, k, 'const');

%% Get element sizes
[c, n, w] = geo2gauss(mesh, 1);
%plot_mesh(mesh);
%plot_elem_normals(mesh);
%%
q = ones(size(Hbm,1),1);

% THIS is OK for alpha = 0
%p_bm = (Hbm - (1/2 + 1i*k/2*alpha)*eye(size(Hbm))) \ ...
%       (Gbm * q + (alpha/2 - 1i/2/k) * q);

p_bm = (Hbm - (1 + alpha*1i*k/2)*eye(size(Hbm))) \ ...
       (Gbm * q - (alpha/2 + 1i/(2*k)) * q);


p = (H - 0.5*eye(size(H)))  \ (G*q);   
pbm = (Hbm - 0.5*eye(size(Hbm))) \ (Gbm * q - 1i/(2*k)*q );
%%
subplot(1,2,1);
plot_mesh(mesh,abs(p_bm)); colorbar;
subplot(1,2,2);
plot_mesh(mesh, abs(p)); colorbar;
