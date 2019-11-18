clear;
close all;

%% parameters
Lx = 1;
Ly = .6;
Le = 2e-2;

%% Build mesh
field = create_slab([Lx, Ly], ceil([Lx Ly] / Le));
bou = drop_mesh_IDs(drop_unused_nodes(flip_mesh(get_boundary(field))));

[xs, ns] = centnorm(bou);
xs = xs(:,1:2);
ns = ns(:,1:2);
xf = centnorm(field);
xf = xf(:,1:2);

%% Create BC
syms x y;
a = 1;
b = 1;
c = 1;
d = 1;
f = a * x +...
    b * y +...
    c * x * y +...
    d * (x^2 - y^2);
gradf = [diff(f, x) diff(f,y)];

us0 = double(subs(f, {x, y}, {xs(:,1), xs(:,2)}));
uf0 = double(subs(f, {x, y}, {xf(:,1), xf(:,2)}));
qs0 = dot(double(subs(gradf, {x, y}, {xs(:,1), xs(:,2)})), ns, 2);

%% Generate BEM matrices
[s_nodes, s_elem] = extract_core_mesh(bou);
[f_nodes, f_elem] = extract_core_mesh(field);
f_elem(:,1) = 22404; % convert to volume quads

%// call C++ code at wave number k = 4
[Ls, Ms, Lf, Mf] = laplace_2d_coll_mex(s_nodes, s_elem, f_nodes, f_elem);

%% Solve Dirichet Problem
qs = Ls \ (Ms * us0);
ufD = Mf * us0 - Lf * qs;
esD = norm(qs - qs0) / norm(qs0);
efD = norm(ufD - uf0) / norm(uf0);

%% Solve Neumann Problem
us = Ms \ (Ls * qs0);
us = us - mean(us) + mean(us0);
ufN = Mf * us - Lf * qs0;
esN = norm(us - us0) / norm(us0);
efN = norm(ufN - uf0) / norm(uf0);

%%
figure;
subplot(2,2,1);
plot_mesh(field, uf0);
colorbar;
title('Analytical');
axis equal tight;
subplot(2,2,2);
plot_mesh(field, ufD);
colorbar;
title(sprintf('Dirichlet, surf err: %.2x, field err: %.2x', esD, efD));
axis equal tight;
subplot(2,2,4);
plot_mesh(field, ufN);
title(sprintf('Neumann, surf err: %.2x, field err: %.2x', esN, efN));
colorbar;
axis equal tight;
