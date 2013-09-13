%// sphere radiator with radius R = 1, divided into 6 elements along the radius.
radiator = quad2tria(create_sphere_boundary(1, 6));
%// CHIEF-points: element centers of cube boundary
chief_points = create_brick_boundary(.4, 5);

%// extract the mesh description matrices
[r_nodes, r_elem] = extract_core_mesh(radiator);
[c_nodes, c_elem] = extract_core_mesh(chief_points);

%// call C++ code at wave number k \approx pi
k = pi + 1.2e-2;
[Ls, Ms, Lc, Mc, Mts, Ns] =...
    helmholtz_bem_3d_fict(r_nodes, r_elem, c_nodes, c_elem, k);

%// define the incident velocity field and analytical solutions
x0 = [0 0 0];
[r_cent, r_norm] = centnorm(radiator);
[ps_ana, qs_ana] = incident('point', x0, r_cent, r_norm, k);

%// acoustic pressure on the surface and in the field points
ps_conv = Ms \ (Ls * qs_ana);               %// conventional
ps_chief = [Ms; Mc] \ ([Ls; Lc] * qs_ana);  %// CHIEF
coup = -1i/k;   %// coupling constant
ps_bm = (Ms + coup * Ns) \ (Ls * qs_ana + coup * Mts * qs_ana);

% // calculate errors
ps_conv_err = abs(ps_conv./ps_ana - 1);     %// conventional
ps_chief_err = abs(ps_chief./ps_ana - 1);   %// CHIEF
ps_bm_err = abs(ps_bm./ps_ana - 1);         %// BM

% // plot results
figure;

plot_mesh(radiator, log10(ps_conv_err));
text(0,1.3,0, 'Conventional', 'HorizontalAlignment', 'center');
plot_mesh(translate_mesh(radiator, [2.3 0 0]), log10(ps_chief_err));
text(2.3,1.3,0, 'CHIEF', 'HorizontalAlignment', 'center');
plot_mesh(translate_mesh(radiator, [2*2.3 0 0]), log10(ps_bm_err));
text(2*2.3,1.3,0, 'Burton-Miller', 'HorizontalAlignment', 'center');
shading flat; caxis(CAx);
view(2);
c = colorbar('SouthOutside');
xlabel(c, 'Log 10 error of solution');

%// display error information
fprintf('Conventional  log10 error: % .3f \n', log10(mean(ps_conv_err)));
fprintf('CHIEF         log10 error: % .3f \n', log10(mean(ps_chief_err)));
fprintf('Burton-Miller log10 error: % .3f \n', log10(mean(ps_bm_err)));
