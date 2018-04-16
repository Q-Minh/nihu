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

kmax = min(mesh_kmax(mesh));
k = .6 * kmax;

[surf_nodes, surf_elements] = extract_core_mesh(mesh);
[field_nodes, field_elements] = extract_core_mesh(field);

[G, H, Ht, D, Gf, Hf] = helmholtz_bem_3d_matrices_gauss(...
    surf_nodes, surf_elements, field_nodes, field_elements, k);

%% check normal derivative matrices
Ht_appr = (Gf - G) / eps;
D_appr = (Hf - H) / eps;
for e = 1 : nElem
    idx = (e-1)*4+(1:4);
    Ht_appr(idx,idx) = Ht(idx,idx);
    D_appr(idx,idx) = D(idx,idx);
end
rel_err_Ht = abs(Ht_appr - Ht) ./ abs(Ht);
rel_err_D = abs(D_appr - D) ./ abs(D);

%% find eleme pairs where relative D error is high

sel = nan(0,2);

for e = 1 : nElem
    eidx = (e-1)*4+(1:4);
    for f = 1 : nElem
    fidx = (f-1)*4+(1:4);
    d = D(eidx,fidx);
    da = D_appr(eidx,fidx);
    if norm(d-da, 'fro')/norm(d, 'fro') > 10
        sel(end+1,1:2) = [e f];
    end
    end
end

for i = 1 : size(sel)
    plot_mesh(mesh, 'elem', sel(i,:)', [1 1]');
    pause;
end

%%
pmesh = mesh;
pmesh.Elements = pmesh.Elements(1:8,:);

figure;
plot_mesh(pmesh, rel_err_Ht(1:8,1));
axis equal;
title(sprintf('Error: %x\n', norm(Ht_appr - Ht, 'fro') / norm(Ht, 'fro')));

figure;
plot_mesh(pmesh, rel_err_D(1:8,1));
axis equal;
title(sprintf('Error: %x\n', norm(D_appr - D, 'fro') / norm(D, 'fro')));

%% Compute pressure
x0 = [1 1 1];
[p0, q0] = incident('point', x0, xc, nc, k); 

N = size(G,2);

I = eye(size(H));
ps = (H - .5 * I) \ (G * q0);

% p0 = -R / (1+1i*k*R) * ones(N,1);
% p0 = ps;

ps_2 = D \ ((Ht + .5 * I) * q0);
alpha = 1i/k;
ps_bm = ((H - .5 * I) + alpha * D) \ ((G + alpha * (Ht + .5 * I)) * q0);
err = norm(ps - p0) / norm(p0);
err_2 = norm(ps_2 - p0) / norm(p0);
err_bm = norm(ps_bm - p0) / norm(p0);

disp([err err_2 err_bm])

%% replace derivative matrices with approximations
ps_2_appr = D_appr \ ((Ht_appr + .5 * I) * q0);
alpha = 1i/k;
ps_bm_appr = ((H - .5 * I) + alpha * D_appr) \ ((G + alpha * (Ht_appr + .5 * I)) * q0);
err_2_appr = norm(ps_2_appr - p0) / norm(p0);
err_bm_appr = norm(ps_bm_appr - p0) / norm(p0);

disp([err err_2_appr err_bm_appr])
