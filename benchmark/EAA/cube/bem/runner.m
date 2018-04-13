clear;

%%
generate_problem;

%% run
system('bem_const.exe data/cube_quad.off data/points.off data/quad_const.xct data/quad_const_ps_0.res data/quad_const_pf_0.res');
system('bem_gauss.exe data/cube_quad.off data/points.off data/quad_gauss.xct data/quad_gauss_ps_0.res data/quad_gauss_pf_0.res');
system('bem_const_bm.exe data/cube_quad.off data/points.off data/quad_const.xct data/quad_const_ps_bm.res data/quad_const_pf_bm.res');
system('bem_gauss_bm.exe data/cube_quad.off data/points.off data/quad_gauss.xct data/quad_gauss_ps_bm.res data/quad_gauss_pf_bm.res');

%% import
load data/freqs

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

data = load('data/quad_gauss_ps_bm.res', '-ascii');
ps_gb = complex(data(:,1:2:end), data(:,2:2:end)).';
data = load('data/quad_gauss_pf_bm.res', '-ascii');
pf_gb = complex(data(:,1:2:end), data(:,2:2:end)).';

%%
figure;
formatfig();
idx = 2;
plot(freqvec, 20*log10(abs(pf_c0(idx,:)/2e-5)), ...
    freqvec, 20*log10(abs(pf_cb(idx,:)/2e-5)), ...
    freqvec, 20*log10(abs(pf_g0(idx,:)/2e-5)), ...
    freqvec, 20*log10(abs(pf_gb(idx,:)/2e-5)));
xlabel('Frequency [Hz]');
ylabel('Sound pressure [dB]');
legend('constant', 'constan bm', 'gauss', 'gauss bm', 'location', 'SouthEast');
setfig('FontSize', 12, 'LineWidth', 1);
title(sprintf('Field Point Pressure in point %d', idx));

%%
figure;
formatfig();
idx = 100;
plot(freqvec, 20*log10(abs(ps_c0(idx,:)/2e-5)), ...
    freqvec, 20*log10(abs(ps_cb(idx,:)/2e-5)), ...
    freqvec, 20*log10(abs(mean(ps_g0(4*idx-3:4*idx,:),1)/2e-5)), ...
    freqvec, 20*log10(abs(mean(ps_gb(4*idx-3:4*idx,:),1)/2e-5)));
xlabel('Frequency [Hz]');
ylabel('Sound pressure [dB]');
legend('constant', 'constan bm', 'gauss', 'gauss bm', 'location', 'SouthEast');
setfig('FontSize', 12, 'LineWidth', 1);
title(sprintf('Surface Pressure in point %d', idx));
