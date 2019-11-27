clear;
close all;

freqvec = .5 : .5 : 300;

p_fmm = nan(length(freqvec), 6);
p_conv = nan(length(freqvec), 6);

le_conv = 100e-3;

for iF = 1 : length(freqvec)
    freq = freqvec(iF);
    
    try
        fname = sprintf('data_conv/gauss_%03dmm/gauss_%03dmm_%gHz_pf.res', ...
            round(le_conv*1000), round(le_conv*1000), freq);
        fid = fopen(fname, 'rt');
        data = fscanf(fid, '%g', 2);
        data = fscanf(fid, '%g', [2 6]);
        p_conv(iF, :) = complex(data(1,:), data(2,:));
        fclose(fid);
    catch
    end
    
    try
        fname = sprintf('data_fmm/conv_data/gauss_quad_bm_10cm_%gpf.res', freq);
        fid = fopen(fname, 'rt');
        data = fscanf(fid, '%g', 1);
        data = fscanf(fid, '%g', [2 6]);
        p_fmm(iF, :) = complex(data(1,:), data(2,:));
        fclose(fid);
    catch
    end
end

%%
close all;
for i = 1 : 6
    figure;
    plot(freqvec, 20*log10(abs(p_fmm(:,i))/2e-5), ...
        freqvec, 20*log10(abs(p_conv(:,i))/2e-5));
    legend('fmm', 'conv');
end
