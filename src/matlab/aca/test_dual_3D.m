%TEST_DUAL_3D test 3D ACA with different meshes

clear;

n = 50;
source = create_circle_boundary(1, n);
source = extrude_mesh(source, [0 0 2*pi/n], n);
receiver = scale_mesh(source, [.5 .5 1]);
xs = centnorm(source);
xr = centnorm(receiver);
nLeaf = 100;
Admit = @(C1, C2)is_admissible_dist(C1, C2, .6);
k = min(mesh_kmax(source));
M = @(i,j)helmholtz_matrix(xr(i,:),xs(j,:),k);
exc = ones(size(xs,1),1);
eps = 1e-2;

[times, result] = aca_tester(M, xr, xs, exc, nLeaf, Admit, eps);

logeps = log10(abs(result.resp./result.resp0-1));

figure;
plot_mesh(receiver, logeps);
c = colorbar;
ylabel(c, 'log10 error [-]');
caxis([-8 0]);
set(gcf, 'renderer', 'painters');
