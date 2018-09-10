clear;

%%
generate_problem;

%% run
% system('bem_const.exe data/sphere_quad.off data/points.off data/quad_const.xct data/quad_const_ps_0.res data/quad_const_pf_0.res');
% system('bem_gauss.exe data/sphere_quad.off data/points.off data/quad_gauss.xct data/quad_gauss_ps_0.res data/quad_gauss_pf_0.res');
system('bem_const_bm.exe data/sphere_quad.off data/points.off data/quad_const.xct data/quad_const_ps_bm.res data/quad_const_pf_bm.res');
system('bem_gauss_bm.exe data/sphere_quad.off data/points.off data/quad_gauss.xct data/quad_gauss_ps_bm.res data/quad_gauss_pf_bm.res');

%% import
R = 3;
r = 4.5;
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

c = 340;
om = 2*pi*freqvec;
k = om/c;
ps_anal = -R ./ (1 + 1i*k*R);
pf_anal = ps_anal .* exp(-1i*k*(r-R))/r*R;

%%
figure;
formatfig();
idx = 2;
plot(k*R, 20*log10(abs(pf_c0(idx,:)/2e-5)), ...
    k*R, 20*log10(abs(pf_cb(idx,:)/2e-5)), ...
    k*R, 20*log10(abs(pf_g0(idx,:)/2e-5)), ...
    k*R, 20*log10(abs(pf_gb(idx,:)/2e-5)), ...
    k*R, 20*log10(abs(pf_anal/2e-5)));
xlabel('Helmholtz number [Hz]');
ylabel('Sound pressure [dB]');
legend('constant', 'constan bm', 'gauss', 'gauss bm', 'anal', 'location', 'SouthEast');
setfig('FontSize', 12, 'LineWidth', 1);
title(sprintf('Field Point Pressure in point %d', idx));

%%
figure;
formatfig();
idx = 2;
plot(k*R, 20*log10(abs(pf_c0(idx,:)./pf_anal)), ...
    k*R, 20*log10(abs(pf_cb(idx,:)./pf_anal)), ...
    k*R, 20*log10(abs(pf_g0(idx,:)./pf_anal)), ...
    k*R, 20*log10(abs(pf_gb(idx,:)./pf_anal)));
xlabel('Helmholtz number [-]');
ylabel('Error [dB]');
legend('constant', 'constan bm', 'gauss', 'gauss bm', 'location', 'SouthEast');
setfig('FontSize', 12, 'LineWidth', 1);
title(sprintf('Field Point Pressure in point %d', idx));
ylim([-1 1]);
grid;

%%
figure;
formatfig();
idx = 2;
plot(freqvec, 20*log10(abs(pf_c0(idx,:)/2e-5)), ...
    freqvec, 20*log10(abs(pf_cb(idx,:)/2e-5)), ...
    freqvec, 20*log10(abs(pf_g0(idx,:)/2e-5)), ...
    freqvec, 20*log10(abs(pf_gb(idx,:)/2e-5)), ...
    freqvec, 20*log10(abs(pf_anal/2e-5)));
xlabel('Frequency [Hz]');
ylabel('Sound pressure [dB]');
legend('constant', 'constan bm', 'gauss', 'gauss bm', 'anal', 'location', 'SouthEast');
setfig('FontSize', 12, 'LineWidth', 1);
title(sprintf('Field Point Pressure in point %d', idx));


%%
figure;
formatfig();
idx = 100;
plot(freqvec, 20*log10(abs(ps_c0(idx,:)./ps_anal)), ...
    freqvec, 20*log10(abs(ps_cb(idx,:)./ps_anal)), ...
    freqvec, 20*log10(abs(mean(ps_g0(4*idx-3:4*idx,:),1)./ps_anal)), ...
    freqvec, 20*log10(abs(mean(ps_gb(4*idx-3:4*idx,:),1)./ps_anal)));
xlabel('Frequency [Hz]');
ylabel('Error [dB]');
legend('constant', 'constan bm', 'gauss', 'gauss bm', 'location', 'SouthEast');
setfig('FontSize', 12, 'LineWidth', 1);
title(sprintf('Surface Pressure in point %d', idx));
ylim([-1 1]);
grid;
