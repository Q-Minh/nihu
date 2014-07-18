clear;

%%
mesh = create_sphere_boundary(1, 6);
m1 = drop_unused_nodes(mesh_section(mesh, [-Inf, -Inf, -Inf; +Inf, 0+1e-3, +Inf]));
m2 = drop_unused_nodes(mesh_section(mesh, [-Inf, 0-1e-3, -Inf; +Inf, +Inf, +Inf]));
m2 = quad2tria(m2);
mesh = merge_coincident_nodes(join_meshes(m1, m2));
k = min(bemkmax(mesh, 8))*.7;
r0 = [.2 .1 0];

%%
tic;
B = ibemB(mesh, k, 'const');
toc;

%%
[xgi, ngi, wgi, igi] = geo2gauss(mesh, 2);
pg = incident('point', r0, xgi, [], k);

%%
tic;
nE = size(mesh.Elements,1);
a = zeros(nE, 1);
for i = 1 : nE
    indi = find(igi == i);
    wi = wgi(indi);
    a(i) = wi .' * pg(indi);
end
toc;

sigma = -B \ a;

%%
field = create_line([2 0 0; 6 0 0], 60);
field = revolve_mesh(field, [0 0 0], [0 0 1], pi/40, 20);

xF = field.Nodes(:,2:end);
pf0 = incident('point', r0, xF, [], k);
nF = size(field.Nodes,1);

tic;
Bf = ibemB(mesh, k, 'lin', xF);
toc;

pf = -Bf * sigma;

%%
figure; plot_mesh(field, imag(pf)); shading interp; colorbar;
figure; plot_mesh(field, imag(pf0)); shading interp; colorbar;
figure; plot_mesh(field, log10(abs(pf./pf0-1))); shading interp; colorbar;