clear;
mesh = create_line(1, 100);
sigma = 2;
d = .5;
[nodes, elements] = extract_core_mesh(mesh);
elements(:,1) = 11202; 	% convert to volume elements
[D, B] = kl1d(nodes, elements, sigma, d);
