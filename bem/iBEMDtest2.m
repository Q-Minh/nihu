clear;

%%
mesh = translate_mesh(create_slab(1, 20), [-.5 -.5]);

mesh = quad2tria(mesh);
% mesh = merge_coincident_nodes(join_meshes(m1, m2));
k = min(bemkmax(mesh, 8))*.7;

%%
% mesh.Elements = mesh.Elements(1,:);
tic;
D = ibemD(mesh, k, 'lin');
toc;

%%
x = mesh.Nodes(:,2);
y = mesh.Nodes(:,3);
nx = 2;
ny = 1;
vs = cos(nx*pi*x) .* cos(ny*pi*y);

%%
M = ibemRHS(mesh);
mu = D \ vs;

%%
field = translate_mesh(rotate_mesh(create_slab(2, 40), [0 0 0], [1 0 0], pi/2), [-1 0 .2]);
field = join_meshes(field, reflect_mesh(field, [0 0 0], [0 0 1]));
xF = field.Nodes(:,2:end);

Df = ibemD(mesh, k, 'lin', xF);

pf = Df * mu;

%%
figure;
plot_mesh(field, real(pf)); shading interp; colorbar;
plot_mesh(mesh);
% figure; plot_mesh(field, log10(abs(pf./pf0-1))); shading interp; colorbar;