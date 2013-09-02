%// sphere radiator with radius R = 1, divided into 8 elements along the radius.
radiator = quad2tria(create_sphere_boundary(1, 6));
chief_points = create_brick_boundary(.4, 5);

%// extract the mesh description matrices
[r_nodes, r_elem] = extract_Boonen_mesh(radiator);
[c_nodes, c_elem] = extract_Boonen_mesh(chief_points);

%// call C++ code at wave number k \approx pi
k = pi + 1.2e-2;
[Ls, Ms, Lc, Mc, Mts, Ns] =...
    helmholtz_bem_3d_fict(r_nodes, r_elem, c_nodes, c_elem, k);

%// define the incident velocity field and analytical solutions
x0 = [-.2 0 .3];
[r_cent, r_norm] = centnorm(radiator);
[ps_ana, qs_ana] = incident('point', x0, r_cent, r_norm, k);

%// acoustic pressure on the surface and in the field points
ps_conv = Ms \ (Ls * qs_ana);
ps_chief = [Ms; Mc] \ ([Ls; Lc] * qs_ana);
coup = -1i/k;
ps_bm = (Ms + coup * Ns) \ (Ls * qs_ana + coup * Mts * qs_ana);

% // calculate errors
ps_conv_err = abs(ps_conv./ps_ana - 1);
ps_chief_err = abs(ps_chief./ps_ana - 1);
ps_bm_err = abs(ps_bm./ps_ana - 1);

% % // plot results
% figure;
% subplot(1,2,1);
% plot_mesh(radiator, real(ps_num));
% plot_mesh(field, real(pf_num)); view(2);  shading flat;
% c1 = colorbar('SouthOutside');
% xlabel(c1, 'Real part of numerical solution');
% subplot(1,2,2);
% plot_mesh(radiator, log10(ps_err));
% plot_mesh(field, log10(pf_err)); view(2); shading flat;
% c2 = colorbar('SouthOutside');
% xlabel(c2, 'Log 10 error of solution');
% 
% %// display error information
% fprintf('Surface log10 error: %f \n', log10(mean(ps_err)));
% fprintf('Field   log10 error: %f \n', log10(mean(pf_err)));
