%TEST_ACA_DUAL_1D test 1D ACA with different meshes

clear;

Ns = 1000;
Nr = 3000;
source = create_line(1, Ns);
receiver = translate_mesh(create_line(1, Nr), [1.5 0 0]);
xs = centnorm(source);
xs = xs(:,1:2);
xr = centnorm(receiver);
xr = xr(:,2:-1:1);

nLeaf = 25;
eta = .8;
Admit = @(C1, C2)is_admissible_dist(C1, C2, eta);
k = min(mesh_kmax(source));
% k = 10;
exc = rand(Ns,1);
eps = 1e-3;

[times, result] = aca_tester_mex(k, xr, xs, exc, nLeaf, Admit, eps);

% figure;
% plot(xr(:,1), log10(abs(result.resp./result.resp0-1)));
% xlabel('x [m]');
% ylabel('log10 error [-]');
% ylim([-10 0]);

times
