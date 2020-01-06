clear;
close all;

freqvec = .5 : .5 : 800;

pf_fmm = nan(length(freqvec), 6);
pf_conv = nan(length(freqvec), 6);

n_fmm = 0;
n_conv = 0;

le_conv = 100e-3;
n_freq = length(freqvec);
%%
for iF = 1 : n_freq
    freq = freqvec(iF);
    
    try
        fname = sprintf('data_conv/gauss_%03dmm/gauss_%03dmm_%gHz_pf.res', ...
            round(le_conv*1000), round(le_conv*1000), freq);
        fid = fopen(fname, 'rt');
        data = fscanf(fid, '%g', 2);
        data = fscanf(fid, '%g', [2 6]);
        pf_conv(iF, :) = complex(data(1,:), data(2,:));
        fclose(fid);
        n_conv = n_conv+1;
    catch
    end
    
    try
        fname = sprintf('data_fmm/gauss_%03dmm/gauss_%03dmm_%gHz_pf.res', ...
            round(le_conv*1000), round(le_conv*1000), freq);
        fid = fopen(fname, 'rt');
        data = fscanf(fid, '%g', 2);
        data = fscanf(fid, '%g', [2 6]);
        pf_fmm(iF, :) = complex(data(1,:), data(2,:));
        fclose(fid);
        n_fmm = n_fmm+1;
    catch
    end
    
    progbar(1, length(freqvec), iF);
end

fprintf('Loaded %4d / %4d CONV field results \n', n_conv, n_freq);
fprintf('Loaded %4d / %4d FMM  field results \n', n_fmm, n_freq);

%%
close all;
for i = 1 : 6
    figure;
    plot(freqvec, 20*log10(abs(pf_fmm(:,i))/2e-5), ...
        freqvec, 20*log10(abs(pf_conv(:,i))/2e-5));
    legend('fmm', 'conv');
end

%% Load surface results
P_conv = nan(length(freqvec), 1);
iters_conv = nan(length(freqvec), 1);

P_fmm = nan(length(freqvec), 1);
iters_fmm = nan(length(freqvec), 1);
n_conv = 0;
n_fmm = 0;
for iF = 1 : length(freqvec) 
    freq = freqvec(iF);
    
    try
        fname = sprintf('data_conv/gauss_%03dmm/gauss_%03dmm_%gHz_ps.res', ...
            round(le_conv*1000), round(le_conv*1000), freq);
        fid = fopen(fname, 'rt');
        data = fscanf(fid, '%g', 2);
        data = fscanf(fid, '%g', [2 data(2)]);
        ps_conv = complex(data(1,:), data(2,:)).';
        data = fscanf(fid, '%g', 1);
        iters_conv(iF) = data(1);
        fclose(fid);
        n_conv = n_conv+1;
        P_conv(iF) = sum(real(ps_conv), 1)*le_conv^2 / 4;
        
        
    catch
    end
    
    try
        fname = sprintf('data_fmm/gauss_%03dmm/gauss_%03dmm_%gHz_ps.res', ...
            round(le_conv*1000), round(le_conv*1000), freq);
        fid = fopen(fname, 'rt');
        data = fscanf(fid, '%g', 2);
        data = fscanf(fid, '%g', [2 data(2)]);
        ps_fmm = complex(data(1,:), data(2,:)).';
        data = fscanf(fid, '%g', 1);
        iters_fmm(iF) = data(1);
        fclose(fid);
        n_fmm = n_fmm+1;
        P_fmm(iF) = sum(real(ps_fmm), 1)*le_conv^2 / 4;
        
        
    catch
    end
    
    %progbar(1, length(freqvec), iF);
    
end
fprintf('Loaded %4d / %4d CONV surface results\n', n_conv, n_freq);
fprintf('Loaded %4d / %4d FMM surface results\n', n_fmm, n_freq);
%%
if (any(P_conv < 0))
    warning('Non-negative radiated sound power found!');
end
%%
figure;
plot(freqvec, P_conv, freqvec, P_fmm);
xlabel('Frequency [Hz]');
ylabel('Radiated sound power');
%%
figure;
plot(freqvec, iters_conv, freqvec, iters_fmm);

%%
mesh_name = sprintf('radiatterer_%03dmm_quad.off', round(1000*le_conv));
mesh = import_off_mesh(fullfile('data', mesh_name));

p_plot_fmm = mean(reshape(ps_fmm, 4, []), 1).';
p_plot_conv = mean(reshape(ps_conv, 4, []), 1).';

figure;
subplot(1,2,1);
plot_mesh(mesh, abs(p_plot_fmm));
axis equal;
subplot(1,2,2);
plot_mesh(mesh, abs(p_plot_conv));
axis equal;