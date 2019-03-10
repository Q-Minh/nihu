rho = 1.2041;
c = 343.21;

Le = 4e-3;
[pac, field, xct] = create_pac_man(Le);

%%
surf_off_name = 'pac_man_surf_4mm.off';
export_off_mesh(pac, surf_off_name);

line_field = field;
line_field.Elements(:,2) = ShapeSet.LinearLine.Id;
line_field.Elements(:,5:6) = field.Elements(:,[5 7]);
line_field.Elements(:,7:end) = [];

field_off_name = 'pac_man_field_4mm.off';
export_off_mesh(line_field, field_off_name);

%%
k = .7 * min(mesh_kmax(pac));
om = k * c;
f = floor(om/2/pi);


pattern = sprintf('pac_man_4mm_%d', f);

%%
export_excitation(xct * -1i*om*rho, k, sprintf('%s.xct', pattern));

command = sprintf('solve_pac_man.exe %s %s %s', surf_off_name, field_off_name, pattern);
[status, result] = system(command);
disp(result);

%%
[ps, ks] = import_response(sprintf('%s_surf.res', pattern));
[pf, kf] = import_response(sprintf('%s_field.res', pattern));

figure;
plot([real(ps) imag(ps)]);

%%
figure;
plot_mesh(field, 20*log10(abs(pf/2e-5)));
shading flat;
axis equal tight off;
cb = colorbar;
ylabel(cb, 'Sound pressure level [dB]');
caxis([120 150]);

figure;
plot_mesh(field, real(pf));
shading flat;
axis equal tight off;
cb = colorbar;
ylabel(cb, 'Sound pressure [Pa]');
hold on;
plot_mesh(pac);
% caxis([120 150]);

