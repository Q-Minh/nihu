clear;
mesh = quad2tria(create_sphere_boundary(1, 10));
var = [2 1; 1 4];
d = .5;
[nodes, elements] = extract_core_mesh(mesh);
[D, B] = kltria_surface(nodes, elements, var, d);
D = (D + D') / 2;
B = (B + B') / 2;
[g, Lambda] = eigs(D, B, 50);
lambda = diag(Lambda);

figure;
n = 20;
subplot(1,2,1);
plot_mesh(mesh, g(1:2:end,n));
subplot(1,2,2);
plot_mesh(mesh, g(2:2:end,n));
