%// sphere radiator with radius R = 1, divided into 6 elements along the radius.
radiator = quad2tria(create_sphere_boundary(1, 6));
field = quad2tria(create_sphere_boundary(2, 6));

%// extract the mesh description matrices
[r_nodes, r_elem] = extract_core_mesh(radiator);
[f_nodes, f_elem] = extract_core_mesh(field);

%// call C++ code
[Ls, Ms, Lf, Mf] = laplace_bem_3d(r_nodes, r_elem, f_nodes, f_elem);

%// define the incident derivative field and analytical solutions
x0 = [0 .2 0];
[r_cent, r_norm] = centnorm(radiator);
[f_cent, f_norm] = centnorm(field);
[ps_ana, qs_ana] = incident('point', x0, r_cent, r_norm, 0);
[pf_ana, qf_ana] = incident('point', x0, f_cent, f_norm, 0);

%// potential on the surface
I = eye(size(Ls));

%// solve Neumann problem
ps_num = (Ms - .5 * I) \ (Ls * qs_ana);
pf_neu_num = Mf * ps_num - Lf * qs_ana;
ps_err = abs(ps_num ./ ps_ana - 1);
pf_neu_err = abs(pf_neu_num ./ pf_ana - 1);
