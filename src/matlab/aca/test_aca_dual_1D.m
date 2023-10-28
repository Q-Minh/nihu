%TEST_ACA_DUAL_1D test 1D ACA with different meshes

clear;

Ns = 100;
Nr = 200;
source = create_line(1, Ns);
receiver = translate_mesh(create_line(1, Nr), [3 0 0]);
xs = centnorm(source);
xs = xs(:,1);
xr = centnorm(receiver);
xr = xr(:,1);

nLeaf = 25;
eta = .8;
Admit = @(C1, C2)is_admissible_dist(C1, C2, eta);
k = min(mesh_kmax(source));
exc = ones(Ns,1);
eps = 1e-2;

[times, result] = aca_tester_mex(k, xr, xs, exc, nLeaf, Admit, eps);

figure;
plot(xr(:,1), log10(abs(result.resp./result.resp0-1)));
xlabel('x [m]');
ylabel('log10 error [-]');
ylim([-10 0]);

times
