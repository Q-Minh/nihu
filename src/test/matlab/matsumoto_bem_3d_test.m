%// sphere radiator with radius R = 1, divided into 6 elements along the radius.
radiator = quad2tria(create_sphere_boundary(1, 4));

%// extract the mesh description matrices
[r_nodes, r_elem] = extract_core_mesh(radiator);

%// call C++ code at wave number k \approx pi
k = pi + 1.2e-2;
[Ls, Ms, Mts, Ns] =...
    matsumoto_bem_3d(r_nodes, r_elem, k);

%// define the incident velocity field and analytical solutions
x0 = [0 0 0];
[r_cent, r_norm] = centnorm(radiator);
[ps_ana, qs_ana] = incident('point', x0, r_cent, r_norm, k);

%// acoustic pressure on the surface and in the field points
ps_conv = Ms \ (Ls * qs_ana);               %// conventional
coup = -1i/k;   %// coupling constant
ps_bm = (Ms + coup * Ns) \ (Ls * qs_ana + coup * Mts * qs_ana);

% // calculate errors
ps_conv_err = abs(ps_conv./ps_ana - 1);     %// conventional
ps_bm_err = abs(ps_bm./ps_ana - 1);         %// BM

% // plot results
figure;

plot_mesh(radiator, log10(ps_conv_err));
text(0,1.3,0, 'Conventional', 'HorizontalAlignment', 'center');
plot_mesh(translate_mesh(radiator, [1*2.3 0 0]), log10(ps_bm_err));
text(1*2.3,1.3,0, 'Burton-Miller', 'HorizontalAlignment', 'center');
shading flat;
view(2);
c = colorbar('SouthOutside');
xlabel(c, 'Log 10 error of solution');

%// display error information
fprintf('Conventional            log10 error: % .3f \n', log10(mean(ps_conv_err)));
fprintf('Burton-Miller Matsumoto log10 error: % .3f \n', log10(mean(ps_bm_err)));
