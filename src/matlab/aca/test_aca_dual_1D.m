%TEST_ACA_DUAL_1D test 1D ACA with different meshes

clear;

Ns = 1000;
Nr = 2000;
source = create_line(1, Ns);
receiver = translate_mesh(create_line(1, Nr), [.5 0 0]);
xs = centnorm(source);
xs = xs(:,1);
xr = centnorm(receiver);
xr = xr(:,1);
nLeaf = 25;
Admit = @(C1, C2)is_admissible_dist(C1, C2, .8);
k = min(mesh_kmax(source));
M = @(i,j)helmholtz_kernel(xr(i,:), xs(j,:), k);
exc = ones(Ns,1);
eps = 1e-2;

[times, result] = aca_tester(M, xr, xs, exc, nLeaf, Admit, eps);

figure;
plot(xr(:,1), log10(abs(result.resp./result.resp0-1)));
xlabel('x [m]');
ylabel('log10 error [-]');
ylim([-10 0]);
