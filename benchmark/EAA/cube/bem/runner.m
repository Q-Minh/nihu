clear;

%%
generate_problem;

%% run
system('bem_const.exe data/cube_quad.off data/points.off data/quad_const.xct data/quad_const_ps_0.res data/quad_const_pf_0.res');
system('bem_gauss.exe data/cube_quad.off data/points.off data/quad_gauss.xct data/quad_gauss_ps_0.res data/quad_gauss_pf_0.res');
system('bem_const_bm.exe data/cube_quad.off data/points.off data/quad_const.xct data/quad_const_ps_bm.res data/quad_const_pf_bm.res');
system('bem_gauss_bm.exe data/cube_quad.off data/points.off data/quad_gauss.xct data/quad_gauss_ps_bm.res data/quad_gauss_pf_bm.res');

%% import
freqs = .5 : .5 : 100;

data = load('data/quad_const_ps_0.res', '-ascii');
ps_c0 = complex(data(:,1:2:end), data(:,2:2:end)).';
data = load('data/quad_const_pf_0.res', '-ascii');
pf_c0 = complex(data(:,1:2:end), data(:,2:2:end)).';

data = load('data/quad_gauss_ps_0.res', '-ascii');
ps_g0 = complex(data(:,1:2:end), data(:,2:2:end)).';
data = load('data/quad_gauss_pf_0.res', '-ascii');
pf_g0 = complex(data(:,1:2:end), data(:,2:2:end)).';

data = load('data/quad_const_ps_bm.res', '-ascii');
ps_cb = complex(data(:,1:2:end), data(:,2:2:end)).';
data = load('data/quad_const_pf_bm.res', '-ascii');
pf_cb = complex(data(:,1:2:end), data(:,2:2:end)).';

%%
figure;
formatfig();
idx = 2;
plot(freqs, 20*log10(abs(pf_cb(idx,:)/2e-5)), ...
    freqs, 20*log10(abs(pf_g0(idx,:)/2e-5)), ...
    freqs, 20*log10(abs(pf_c0(idx,:)/2e-5)));
xlabel('Frequency [Hz]');
ylabel('Sound pressure [dB]');
legend('constant burton', 'gauss', 'constant', 'location', 'SouthEast');
setfig('FontSize', 12, 'LineWidth', 1);
title(sprintf('Field Point Pressure in point %d', idx));

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
title(sprintf('Surface Pressure in point %d', idx));

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