clear
radiator = quad2tria(create_sphere_boundary(1, 5));
field_mesh = create_sphere_boundary(2, 5);
fp = centnorm(field_mesh);

[r_nodes, r_elements] = extract_Boonen_mesh(radiator);
[f_nodes, f_elements] = extract_Boonen_mesh(field_mesh);

[r_nodes_old, r_elements_old] = extract_bem_mesh(radiator);
[f_nodes_old, f_elements_old] = extract_bem_mesh(field_mesh);

k = 1;

tic;
[Ls, Ms, Lf, Mf] = acoustic_bem(r_nodes, r_elements, f_nodes, f_elements, k);
t_Booni = toc;

tic;
[Hs, Gs] = bemHG(radiator, k, 'const');
t_old_rad = toc;
tic;
[Hf, Gf] = bemHG(radiator, k, 'const', fp);
t_old_field = toc;

x0 = [.2 .3 0];
[cs, n] = centnorm(radiator);
[ps_ana, qs_ana] = incident('point', x0, cs, n, k);
[cf, n] = centnorm(field_mesh);
pf_ana = incident('point', x0, cf, [], k);

ps_Boonen = Ms \ (Ls * qs_ana);
pf_Boonen = Mf * ps_Boonen - Lf * qs_ana;
err_s_Boonen = abs(ps_Boonen ./ ps_ana - 1);
err_f_Boonen = abs(pf_Boonen ./ pf_ana - 1);

ps_oldSchool = (Hs - .5*eye(size(Hs))) \ (Gs * qs_ana);
pf_oldSchool = Hf * ps_oldSchool - Gf * qs_ana;
err_s_oldSchool = abs(ps_oldSchool ./ ps_ana - 1);
err_f_oldSchool = abs(pf_oldSchool ./ pf_ana - 1);

fprintf(1, 'log10 mean eps | surface, field\n');
fprintf(1, 'Boonen:          %.2f\t%.2f\n', log10(mean(err_s_Boonen)), log10(mean(err_f_Boonen)));
fprintf(1, 'oldSchool:       %.2f\t%.2f\n', log10(mean(err_s_oldSchool)), log10(mean(err_f_oldSchool)));
fprintf(1, 'time gain:\t%.2f\n', t_Booni/(t_old_field+t_old_rad));