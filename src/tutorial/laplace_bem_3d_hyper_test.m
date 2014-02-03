%// sphere radiator with radius R = 1, divided into 6 elements along the radius.
radiator = create_sphere_boundary(1, 6);

%// extract the mesh description matrices
[r_nodes, r_elem] = extract_core_mesh(radiator);

%// call C++ code at wave number k \approx pi
[Ls, Ms, Mts, Ns] = laplace_bem_3d_hyper(r_nodes, r_elem);

%// define the incident velocity field and analytical solutions
x0 = [0 0 0];
[r_cent, r_norm] = centnorm(radiator);
[ps_ana, qs_ana] = incident('point', x0, r_cent, r_norm, 0);

%// acoustic pressure on the surface and in the field points
ps_conv = Ms \ (Ls * qs_ana);       %// conventional
ps_hyper = (Ns) \ (Mts * qs_ana);      %// hypersingular

% // calculate errors
ps_conv_err = abs(ps_conv./ps_ana - 1);     %// conventional
ps_hyper_err = abs(ps_hyper./ps_ana - 1);   %// hypersingular

% // plot results
figure;

plot_mesh(radiator, log10(ps_conv_err));
text(0,1.3,0, 'Conventional', 'HorizontalAlignment', 'center');
plot_mesh(translate_mesh(radiator, [2.3 0 0]), log10(ps_hyper_err));
text(2*2.3,1.3,0, 'Hypersingular', 'HorizontalAlignment', 'center');
shading flat;
view(2);
c = colorbar('SouthOutside');
xlabel(c, 'Log 10 error of solution');

%// display error information
fprintf('Conventional  log10 error: % .3f\n', log10(mean(ps_conv_err)));
fprintf('Hypersingular log10 error: % .3f\n', log10(mean(ps_hyper_err)));
