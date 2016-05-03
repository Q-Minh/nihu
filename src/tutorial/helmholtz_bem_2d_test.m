[pot_bb, pot_path] = read_epspath('potato.eps');
siz = norm(diff(pot_bb));
radiator = meshpath(pot_path, siz/100);
radiator = translate_mesh(scale_mesh(radiator, 1/siz), [.2 .2]);

% radiator = translate_mesh(create_slab_boundary([1 1], [500 500]), [-.5 -.5]);
field = create_slab([1.5 1], [150 100]);

%// extract the mesh description matrices
[r_nodes, r_elem] = extract_core_mesh(radiator);
[f_nodes, f_elem] = extract_core_mesh(field);
%// convert slabs to 2D volume elements
f_elem(:,1) = 22404;

%// call C++ code at wave number k = 4
k = .7 * min(min(mesh_kmax(radiator)), min(mesh_kmax(field)));
[Ls, Ms, Lf, Mf] = helmholtz_bem_2d(r_nodes, r_elem, f_nodes, f_elem, k);

%// define the incident velocity field and analytical solutions
x0 = [0. 0. 0.];
[r_cent, r_norm] = centnorm(radiator);
[ps_ana, qs_ana] = incident('line', x0, r_cent, r_norm, k);

%// analytical solution on the field
f_cent = centnorm(field);
pf_ana = incident('point', x0, f_cent, [], k);

%// acoustic pressure on the surface and in the field points
% ps_num = Ms \ (Ls * qs_ana);
pf_num = Mf*ps_ana - Lf*qs_ana;

% // calculate errors
pf_err = abs(pf_num-pf_ana)./abs(pf_ana);

% // plot results
figure;
hold on;
plot(x0(:,1), x0(:,2), 'k*');
plot_mesh(radiator);
plot_mesh(field, real(pf_num)); view(2);  shading flat;
axis equal tight;
c1 = colorbar('SouthOutside');
xlabel(c1, 'Real part of numerical solution');

