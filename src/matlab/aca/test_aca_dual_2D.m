%TEST_ACA_DUAL_2D test 2D ACA with different meshes

clear;

source = create_circle_boundary(1, 1000);
receiver = translate_mesh(source, [1 0 0]);
xs = centnorm(source);
xs = xs(:,1:2);
xr = centnorm(receiver);
xr = xr(:,1:2);
nLeaf = 25;
Admit = @(C1, C2)is_admissible_dist(C1, C2, .5);
k = min(mesh_kmax(source));
M = @(i,j)helmholtz_kernel(xr(i,:),xs(j,:),k);

exc = ones(size(xs,1),1);
eps = 1e-2;

[times, result] = aca_tester(M, xr, xs, exc, nLeaf, Admit, eps);

logeps = log10(abs(result.resp./result.resp0-1));

figure;
phi = atan2(xr(:,2), xr(:,1));
plot(phi, logeps);
ylabel('log10 error [-]');
ylim([-8 0]);
