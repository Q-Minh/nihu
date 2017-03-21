clear;
mesh = quad2tria(create_sphere_boundary(1, 10));
sigma = 2;
d = .5;
[nodes, elements] = extract_core_mesh(mesh);
[D, B] = kltria_surface(nodes, elements, sigma, d);
[g, Lambda] = eigs(D, B, 50);
lambda = diag(Lambda);