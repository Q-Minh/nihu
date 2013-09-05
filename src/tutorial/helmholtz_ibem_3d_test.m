%// sphere radiator with radius R = 1, divided into 6 elements along the radius.
radiator = quad2tria(create_sphere_boundary(1, 5));
%// field point mesh, a line revolved around the z axis
field = revolve_mesh(create_line([1.125, 0, 0; 3.625, 0, 0], 40), pi/100, 50, [0 0 1]);

%// extract the mesh description matrices
[r_nodes, r_elem] = extract_Boonen_mesh(radiator);
[f_nodes, f_elem] = extract_Boonen_mesh(field);

%// call C++ code at wave number k = 0
k = 4;
[Ms, Mf] = helmholtz_bem_indirect_dirichlet(r_nodes, r_elem, f_nodes, f_elem, k);

%// define the incident velocity field and analytical solutions
x0 = [.4 .3 0];
[r_cent, r_norm] = centnorm(radiator);
[ps_ana, qs_ana] = incident('point', x0, r_cent, r_norm, k);

%// analytical solution on the field
f_cent = centnorm(field);
pf_ana = incident('point', x0, f_cent, [], k);

%// acoustic pressure on the surface and in the field points
sigma_inf = Ms \ ps_ana;
pf_num = Mf * sigma_inf;

% // calculate errors
pf_err = abs(pf_num./pf_ana - 1);

% // plot results
figure;
subplot(1,2,1);
plot_mesh(field, real(pf_num)); view(2);  shading flat;
c1 = colorbar('SouthOutside');
xlabel(c1, 'Real part of numerical solution');
subplot(1,2,2);
plot_mesh(field, log10(pf_err)); view(2);  shading flat;
c1 = colorbar('SouthOutside');
xlabel(c1, 'Log10 error');

%// display error information
fprintf('Field   log10 error: %f \n', log10(mean(pf_err)));
