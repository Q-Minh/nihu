clear;

%%
Le = 1e-2;

[pac, field, k, q_surf, qs_scat_line, pf_in_line] = create_pac_man(Le);

%% export surface mesh
surf_off_name = sprintf('pac_man_surf_%gmm.off', Le*1000);
export_off_mesh(pac, surf_off_name);

%% export field point mesh
line_field = field;
line_field.Elements(:,2) = ShapeSet.LinearLine.Id;
line_field.Elements(:,5:6) = field.Elements(:,[5 7]);
line_field.Elements(:,7:end) = [];

field_off_name = sprintf('pac_man_field_%gmm.off', Le*1000);
export_off_mesh(line_field, field_off_name);

%%
pattern = sprintf('pac_man_%gmm', 1000*Le);

%% solve radiation problem
export_excitation(q_surf, k, sprintf('%s_rad.xct', pattern));
command = sprintf('solve_pac_man.exe %s %s %s', surf_off_name, field_off_name, sprintf('%s_rad', pattern));
[status_rad, result_rad] = system(command);
disp(result_rad);
ps_rad = import_response(sprintf('%s_rad_surf.res', pattern));
pf_rad = import_response(sprintf('%s_rad_field.res', pattern));

figure;
plot([real(ps_rad) imag(ps_rad)]);

%% solve line scattering problem
export_excitation(qs_scat_line, k, sprintf('%s_line.xct', pattern));
command = sprintf('solve_pac_man.exe %s %s %s', surf_off_name, field_off_name, sprintf('%s_line', pattern));
[status_line, result_line] = system(command);
disp(result_line);
ps_line = import_response(sprintf('%s_line_surf.res', pattern));
pf_line = import_response(sprintf('%s_line_field.res', pattern));

figure;
plot([real(ps_line) imag(ps_line)]);

%%
figure;
plot_mesh(field, 20*log10(abs(pf_in_line+pf_line)/2e-5));
shading flat;
axis equal tight off;
cb = colorbar;
ylabel(cb, 'Sound pressure level [dB]');
% caxis([120 150]);

figure;
plot_mesh(field, real(pf_in_line+pf_line));
shading flat;
axis equal tight off;
cb = colorbar;
ylabel(cb, 'Sound pressure [Pa]');
hold on;
plot_mesh(pac);
% caxis([120 150]);

