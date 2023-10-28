clear;
mesh = create_circle_boundary(1, 200);
var = 2;
d = .5;
[nodes, elements] = extract_core_mesh(mesh);
[D, B] = klline_surface(nodes, elements, var, d);
[g, Lambda] = eigs(D, B, 50);
lambda = diag(Lambda);
