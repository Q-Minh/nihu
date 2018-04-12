clear;

%%
generate_problem;

%% run
system('bem_const.exe data/cube_quad.off data/points.off data/quad_const');
system('bem_gauss.exe data/cube_quad.off data/points.off data/quad_gauss');
system('bem_const_bm.exe data/cube_quad.off data/points.off data/quad_const_bm');
system('bem_gauss_bm.exe data/cube_quad.off data/points.off data/quad_gauss_bm');

%% import
[freqs_cb, pf_cb, ps_cb, iters_cb] = import_data('quad_const_bm');
[freqs_gb, pf_gb, ps_gb, iters_gb] = import_data('quad_gauss_bm');
[freqs_c0, pf_c0, ps_c0, iters_c0] = import_data('quad_const_0');
[freqs_g0, pf_g0, ps_g0, iters_g0] = import_data('quad_gauss_0');
ps_gb = squeeze(mean(reshape(ps_gb, 4, [], size(ps_gb,2)), 1));
ps_g0 = squeeze(mean(reshape(ps_g0, 4, [], size(ps_gb,2)), 1));

%%
figure;
formatfig();
idx = 2;
plot(freqs_cb, 20*log10(abs(pf_cb(idx,:)/2e-5)), ...
    freqs_gb, 20*log10(abs(pf_gb(idx,:)/2e-5)), ...
    freqs_c0, 20*log10(abs(pf_c0(idx,:)/2e-5)), ...
    freqs_g0, 20*log10(abs(pf_g0(idx,:)/2e-5)));
xlabel('Frequency [Hz]');
ylabel('Sound pressure [dB]');
legend('constant burton', 'linear burton', 'constant', 'linear', 'location', 'SouthEast');
setfig('FontSize', 12, 'LineWidth', 1);

%%
figure;
formatfig();
idx = 24;
plot(freqs_cb, 20*log10(abs(ps_cb(idx,:)/2e-5)), ...
    freqs_gb, 20*log10(abs(ps_gb(idx,:)/2e-5)), ...
    freqs_c0, 20*log10(abs(ps_c0(idx,:)/2e-5)), ...
    freqs_g0, 20*log10(abs(ps_g0(idx,:)/2e-5)));
xlabel('Frequency [Hz]');
ylabel('Sound pressure [dB]');
legend('constant burton', 'linear burton', 'constant', 'linear', 'location', 'SouthEast');
setfig('FontSize', 12, 'LineWidth', 1);

%%
figure;
formatfig();
idx = 1;
plot(freqs_cb, iters_cb, ...
    freqs_gb, iters_gb, ...
    freqs_c0, iters_c0, ...
    freqs_g0, iters_g0);
xlabel('Frequency [Hz]');
ylabel('#iterations [-]');
legend('constant burton', 'linear burton', 'constant', 'linear', 'location', 'SouthEast');
setfig('FontSize', 12, 'LineWidth', 1);

%%
mesh = import_off_mesh('data/cube_quad.off');
points = import_off_mesh('data/points.off');
figure;
plot_mesh(mesh, real(ps_c0(:,42)));
axis equal tight;
shading flat
hold on;
plot3(points.Nodes(:,2), points.Nodes(:,3), points.Nodes(:,4), 'r*');