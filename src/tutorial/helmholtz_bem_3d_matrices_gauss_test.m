clear;

%% Create surface mesh

R = 1;
% M = 5;
% L = [2 2 2];
% M = [6 6 6];
Le = 1e-1;
mesh = create_radiatterer(Le);
% mesh = create_brick_boundary(L, M);
nElem = size(mesh.Elements,1);

%% Create field point mesh
[xc, nc] = geo2gauss(mesh, 2);

eps = 1e-3;

field = create_empty_mesh();
Nf = size(xc,1);
field.Nodes = [(1 : Nf)'  xc + eps * nc];
field.Elements(1:Nf,1) = 1:Nf;
field.Elements(:,2) = ShapeSet.LinearQuad.Id;
field.Elements(:,3) = 1;
field.Elements(:,4) = 1;
field.Elements(:,5:8) = repmat(field.Nodes(:,1), 1, 4);

f = 20;
c = 340;
k = 2*pi*f/c;

kmax = min(mesh_kmax(mesh));
% k = .6 * kmax;

%% Compute BEM matrices
[surf_nodes, surf_elements] = extract_core_mesh(mesh);
[field_nodes, field_elements] = extract_core_mesh(field);

[G, H, Ht, D, Gf, Hf] = helmholtz_bem_3d_matrices_gauss(...
    surf_nodes, surf_elements, field_nodes, field_elements, k);

%% Compute pressure
% x0 = [1 1 1];
% [p0, q0] = incident('point', x0, xc, nc, k); 
q0 = ones(4*nElem,1);

I = eye(size(H));

%% Solve regular
rhs = G * q0;
A = H - .5 * I;
ps = gmres(A, rhs, 100, 1e-4, 100);

figure;
plot_mesh(mesh, abs(mean(reshape(ps, 4, []), 1)'));
shading flat;
colorbar;
axis equal;
title('Regular');

%% Solve BM
alpha = 1i/k;
rhs = ((G + alpha * (Ht + .5 * I)) * q0);
A = ((H - .5 * I) + alpha * D);
ps_bm = gmres(A, rhs, 100, 1e-4, 100);

figure;
plot_mesh(mesh, abs(mean(reshape(ps_bm, 4, []), 1)'));
shading flat;
colorbar;
axis equal;
title('Burton-Miller');

%% Solve approx BM
Ht_appr = (Gf - G) / eps;
D_appr = (Hf - H) / eps;
for e = 1 : nElem
    idx = (e-1)*4+(1:4);
    Ht_appr(idx,idx) = Ht(idx,idx);
    D_appr(idx,idx) = D(idx,idx);
end

rhs = ((G + alpha * (Ht_appr + .5 * I)) * q0);
A = ((H - .5 * I) + alpha * D_appr);
ps_bm_appr = gmres(A, rhs, 100, 1e-4, 100);

figure;
plot_mesh(mesh, abs(mean(reshape(ps_bm_appr, 4, []), 1)'));
shading flat;
colorbar;
axis equal;
title('Burton-Miller approx');

%% Solve derivative
rhs = (Ht + .5 * I) * q0;
A = D;
ps_2 = gmres(A, rhs, 100, 1e-4, 100);

figure;
plot_mesh(mesh, abs(mean(reshape(ps_2, 4, []), 1)'));
shading flat;
colorbar;
axis equal;
title('Derivative');


