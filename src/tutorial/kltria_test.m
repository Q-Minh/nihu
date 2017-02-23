clear;
mesh = quad2tria(create_slab([.5 1], [20 40]));
sigma = 2;
d = .5;
[nodes, elements] = extract_core_mesh(mesh);
elements(:,1) = 22303; 	% convert to volume triangles
[D, B] = kltria(nodes, elements, sigma, d);
[g, Lambda] = eigs(D, B, 50);
lambda = diag(Lambda);