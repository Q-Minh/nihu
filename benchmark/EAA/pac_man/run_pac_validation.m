clear;

%%
Le = 1e-3;
freq = 2000;
[pac, field, k, q_surf] = create_pac_man_validation(Le, freq);
xf = centnorm(field);
phi = atan2(xf(:,2), xf(:,1));
phi(phi < 0) = phi(phi < 0) + 2*pi;

%% export surface mesh
surf_off_name = sprintf('data/pac_man_surf_%gmm.off', Le*1000);
export_off_mesh(pac, surf_off_name);

%% export field point mesh
field_off_name = sprintf('data/pac_man_field_%gmm.off', Le*1000);
export_off_mesh(field, field_off_name);

%%
pattern = sprintf('data/pac_man_%gmm', 1000*Le);
exe_name = 'helmholtz_2d_wb_fmm_standalone.exe';

%% solve radiation problem
export_excitation(q_surf, k, sprintf('%s_rad.xct', pattern));
command = sprintf('%s %s %s %s', exe_name, surf_off_name, field_off_name, sprintf('%s_rad', pattern));
[status_rad, result_rad] = system(command);
disp(result_rad);
ps_rad = import_response(sprintf('%s_rad_surf.res', pattern));
pf_rad = import_response(sprintf('%s_rad_field.res', pattern));
save(sprintf('%s_rad_surf', pattern), 'result_rad', 'ps_rad', 'pf_rad', 'pac', 'field');

%% Import Excel
csv_name = 'ref_impl\reference_solution_surface_vibration.csv';
data = csvread(csv_name, 1, 1);
% phi_anal = data(:,1);
phi_anal = linspace(0, 360, 72);
pf_rad_anal = data(:,end-1);

%%
figure;
plot([real(ps_rad) imag(ps_rad)]);
title('Surface pressure');

figure;
plot(phi/pi*180, real(pf_rad),...
    phi_anal, real(pf_rad_anal), ...
    phi/pi*180, imag(pf_rad),...
    phi_anal, imag(pf_rad_anal));
legend('real fmm', 'real anal', 'imag fmm', 'imag anal');
title('Field pressure');

