%// circular radiator with radius R = 1, divided into 15 elements along the radius.
radiator = create_circle(1, 15);
%// rectangular field point mesh, defined by 4 corners and division parameters.
field = create_slab([-2 0 .1; 2 0 .1; 2 0 4.1; -2 0 4.1], [40 40]);

%// extract the mesh description matrices
[r_nodes, r_elem] = extract_Boonen_mesh(radiator);
[f_nodes, f_elem] = extract_Boonen_mesh(field);

%// call C++ code at wave number k = 10
[Z_trans, Z_rad] = rayleigh_integral_3d(r_nodes, r_elem, f_nodes, f_elem, 10);

%// constant velocity over the radiator
vn = ones(size(Z_trans,2),1);
%// radiated pressure computed by matrix multiplication
pf = Z_trans * vn;
pr = Z_rad * vn;

%// plot results
figure;
plot_mesh(radiator, 20*log10(abs(pr)));
plot_mesh(field, 20*log10(abs(pf)));

