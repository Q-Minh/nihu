clear;
close all;

c = 340;
n_field = 6;                    % Number of field points

% Cases
cases = { ...
    'conv', 'gauss', 100e-3;
    'fmm',  'gauss', 100e-3;
    'fmm',  'gauss',  50e-3;
    'fmm',  'gauss',  33e-3;
};

n_cases = size(cases, 1);

f = cell(n_cases, 1);       % Frequencies in Hz
pf = cell(n_cases, 1);      % Field point pressures


% Load conventional data
for i_case = 1 : n_cases
    method = cases{i_case, 1};
    form   = cases{i_case, 2};      % Formalism (const / gauss)
    le     = cases{i_case, 3};
    fname_pat = sprintf('data_%s/%s_%03dmm/%s_%03dmm_*Hz_pf.res', ...
       method, form, round(le*1000), form, round(le*1000));
    files = dir(fname_pat);    
    n_files = numel(files);
    
    f{i_case} = nan(n_files, 1);
    pf{i_case} = nan(n_files, n_field);
    
    for i_file = 1 : n_files
        progbar(1, n_files, i_file);
        
        fname = sprintf('data_%s/%s_%03dmm/%s', ...
            method, form, round(le*1000), files(i_file).name);
        fid = fopen(fname, 'rt');
        % Read first two values and drop 2nd
        data = fscanf(fid, '%g', 2);
        % Store freq 
        f{i_case}(i_file) = data(1) * c / (2*pi);
        % Read the field point data
        data = fscanf(fid, '%g', [2 6]);
        % Store complex field point pressures
        pf{i_case}(i_file, :) = complex(data(1,:), data(2,:));
        
        fclose(fid);
    end
    
    % Sort the results
    [f{i_case}, idx] = sort(f{i_case}, 'ascend');
    pf{i_case} = pf{i_case}(idx, :);
    
    fprintf('Loaded %d field point results of case: %s %s le = %d mm\n', ...
        n_files, method, form, round(le*1e3));
end

%%
close all;
for i = 1 : n_field
    figure;
    hold on;
    leg = cell(n_cases, 1);
    for i_case = 1 : n_cases
        plot(f{i_case}, 20*log10(abs(pf{i_case}(:,i))/2e-5));
        leg{i_case} = sprintf('%s %s le = %3d', ...
            cases{i_case, 1}, cases{i_case, 2}, round(cases{i_case, 3}*1e3));
    end
    legend(leg);
end

%% Load surface results
P_conv = nan(length(freqvec), 1);
iters_conv = nan(length(freqvec), 1);

P_fmm = nan(length(freqvec), 1);
iters_fmm = nan(length(freqvec), 1);
n_cases = 0;
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
        n_cases = n_cases+1;
        P_conv(iF) = sum(real(ps_conv), 1)*le_conv^2 / 4;
        
        
    catch
    end
    
    try
        fname = sprintf('data_fmm/gauss_%03dmm/gauss_%03dmm_%gHz_ps.res', ...
            round(le_fmm*1000), round(le_fmm*1000), freq);
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
fprintf('Loaded %4d / %4d CONV surface results\n', n_cases, n_freq);
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

freq = 900;
    
fname = sprintf('data_conv/gauss_%03dmm/gauss_%03dmm_%gHz_ps.res', ...
    round(le_conv*1000), round(le_conv*1000), freq);
fid = fopen(fname, 'rt');
data = fscanf(fid, '%g', 2);
data = fscanf(fid, '%g', [2 data(2)]);
ps_conv = complex(data(1,:), data(2,:)).';

fname = sprintf('data_fmm/gauss_%03dmm/gauss_%03dmm_%gHz_ps.res', ...
    round(le_fmm*1000), round(le_fmm*1000), freq);
fid = fopen(fname, 'rt');
data = fscanf(fid, '%g', 2);
data = fscanf(fid, '%g', [2 data(2)]);
ps_fmm = complex(data(1,:), data(2,:)).';


mesh_name_conv = sprintf('radiatterer_%03dmm_quad.off', round(1000*le_conv));
mesh_conv = import_off_mesh(fullfile('data', mesh_name_conv));

mesh_name_fmm = sprintf('radiatterer_%03dmm_quad.off', round(1000*le_fmm));
mesh_fmm = import_off_mesh(fullfile('data', mesh_name_fmm));

p_plot_fmm = mean(reshape(ps_fmm, 4, []), 1).';
p_plot_conv = mean(reshape(ps_conv, 4, []), 1).';

figure;
subplot(1,2,1);
plot_mesh(mesh_fmm, abs(p_plot_fmm));
axis equal;
subplot(1,2,2);
plot_mesh(mesh_conv, abs(p_plot_conv));
axis equal;