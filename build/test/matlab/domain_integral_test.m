clear all mex classes;

domain = quad2tria(create_slab(1, 10));
boundary = get_boundary(domain);

[nodes, delements] = extract_core_mesh(domain);
delements(:,1) = 20000 + mod(delements(:,1), 10000);
[~, belements] = extract_core_mesh(boundary);

D = domain_integral(nodes, delements, belements);
