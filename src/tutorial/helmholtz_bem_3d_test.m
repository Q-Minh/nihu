%// sphere radiator with radius R = 1, divided into 8 elements along the radius.
radiator = create_sphere_boundary(1, 8);

%// the radiator is divided into two parts to get tria and quad elements
rad_left  = drop_unused_nodes(mesh_section(radiator, [-inf, -inf, -inf; 1e-3, inf, inf]));
rad_right = drop_unused_nodes(mesh_section(radiator, [-1e-3, -inf, -inf; inf, inf, inf]));
radiator = merge_coincident_nodes(join_meshes(rad_left, quad2tria(rad_right)));

%// field point mesh, a line revolved around the z axis
field = revolve_mesh(create_line([1.125, 0, 0; 3.625, 0, 0], 40), pi/100, 50, [0 0 1]);

%// extract the mesh description matrices
[r_nodes, r_elem] = extract_Boonen_mesh(radiator);
[f_nodes, f_elem] = extract_Boonen_mesh(field);

%// call C++ code at wave number k = 4
k = 4;
[Ls, Ms, Lf, Mf] = helmholtz_bem_3d(r_nodes, r_elem, f_nodes, f_elem, k);

%// define the incident velocity field and analytical solutions
x0 = [-.2 0 .3];
[r_cent, r_norm] = centnorm(radiator);
[ps_ana, qs_ana] = incident('point', x0, r_cent, r_norm, k);

%// analytical solution on the field
f_cent = centnorm(field);
pf_ana = incident('point', x0, f_cent, [], k);

%// acoustic pressure on the surface and in the field points
ps_num = Ms \ (Ls * qs_ana);
pf_num = Mf*ps_num - Lf*qs_ana;

% // calculate errors
ps_err = abs(ps_num-ps_ana)./abs(ps_ana);
pf_err = abs(pf_num-pf_ana)./abs(pf_ana);

% // plot results
figure;
subplot(1,2,1);
plot_mesh(radiator, real(ps_num));
plot_mesh(field, real(pf_num)); view(2);  shading flat;
c1 = colorbar('SouthOutside');
xlabel(c1, 'Real part of numerical solution');
subplot(1,2,2);
plot_mesh(radiator, log10(ps_err));
plot_mesh(field, log10(pf_err)); view(2); shading flat;
c2 = colorbar('SouthOutside');
xlabel(c2, 'Log 10 error of solution');

%// display error information
fprintf('Surface log10 error: %f \n', log10(mean(ps_err)));
fprintf('Field   log10 error: %f \n', log10(mean(pf_err)));
