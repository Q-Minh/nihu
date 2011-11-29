clear;

%%
mesh = create_sphere_boundary(1, 6);
% m1 = drop_unused_nodes(mesh_section(mesh, [-Inf, -Inf, -Inf; +Inf, 0+1e-3, +Inf]));
% m2 = drop_unused_nodes(mesh_section(mesh, [-Inf, 0-1e-3, -Inf; +Inf, +Inf, +Inf]));
mesh = quad2tria(mesh);
% mesh = merge_coincident_nodes(join_meshes(m1, m2));
k = min(bemkmax(mesh, 8))*.7;
r0 = [-.7 .2 .5];

%%
% mesh.Elements = mesh.Elements(1,:);
tic;
D = ibemD(mesh, k, 'lin');
toc;

%%
[pg, qg] = incident('point', r0, mesh.Nodes(:,2:4), mesh.Nodes(:,2:4), k);

%%
M = ibemRHS(mesh);
mu = D \ (M * qg);

%%
field = create_line([2 0 0; 6 0 0], 60);
field = revolve_mesh(field, [0 0 0], [0 0 1], pi/40, 20);

xF = field.Nodes(:,2:end);
pf0 = incident('point', r0, xF, [], k);
nF = size(field.Nodes,1);

tic;
Df = ibemD(mesh, k, 'lin', xF);
toc;

pf = Df * mu;

%%
figure; plot_mesh(field, imag(pf)); shading interp; colorbar;
figure; plot_mesh(field, imag(pf0)); shading interp; colorbar;
figure; plot_mesh(field, log10(abs(pf./pf0-1))); shading interp; colorbar;