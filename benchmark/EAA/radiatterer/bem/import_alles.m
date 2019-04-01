clear;

load ../fem/data/result;

pat_fmm = 'data/gauss_bm_quad_05cm_';
pat_conv = 'data/gauss_quad_bm_10cm_';

freqvec = .5 : .5 : 1500;
nF = length(freqvec);

pf_fem = nan(nF, 6);
pf_fem(1:2000, :) = p_field.';
pf_conv = nan(nF, 6);
pf_fmm = nan(nF, 6);

P_fmm = nan(nF,1);
iters_fmm = nan(nF,1);
P_conv = nan(nF,1);
iters_conv = nan(nF,1);


for iF = 1 : length(freqvec)
    progbar(1, length(freqvec), iF);
    freq = freqvec(iF);
    
    fid = fopen(sprintf('%s%gpf.res', pat_fmm, freq), 'rt');
    if fid ~= -1
        fscanf(fid, '%g', [2 1]);
        data = fscanf(fid, '%g', [2 6]);
        pf_fmm(iF,:) = complex(data(1,:), data(2,:));
        fclose(fid);
    end
    
    fid = fopen(sprintf('%s%gps.res', pat_fmm, freq), 'rt');
    if fid ~= -1
        fscanf(fid, '%g', 1);
        n = fscanf(fid, '%g', 1);
        data = fscanf(fid, '%g', [2 n]);
        ps = complex(data(1,:), data(2,:));
        P_fmm(iF) = sum(real(ps));
        iters_fmm(iF) = fscanf(fid, '%g', 1);
        fclose(fid);
    end
    
    fid = fopen(sprintf('%s%gpf.res', pat_conv, freq), 'rt');
    if fid ~= -1
        frq = fscanf(fid, '%g', 1);
        data = fscanf(fid, '%g', [2 6]);
        pf_conv(iF,:) = complex(data(1,:), data(2,:));
        fclose(fid);
    end
    
    
    fid = fopen(sprintf('%s%gps.res', pat_conv, freq), 'rt');
    if fid ~= -1
        fscanf(fid, '%g', 1);
        iters_conv(iF) = fscanf(fid, '%g', 1);
        data = fscanf(fid, '%g', [2 inf]);
        ps = complex(data(1,:), data(2,:));
        P_conv(iF) = sum(real(ps));
        fclose(fid);
    end
end

%%
figure;
formatfig([14, 4], [0 .5 .5 .5]);
n = 4;
plot(freqvec, 20*log10(abs(pf_conv(:,n))/2e-5), ...
    freqvec, 20*log10(abs(pf_fmm(:,n))/2e-5), ...
    freqvec, 20*log10(abs(pf_fem(:,n))/2e-5));
legend({'conv BEM', 'FMM', 'FEM + IEM'}, 'Location', 'NorthWest');
xlabel('Frequency [Hz]');
ylabel('SPL [dB]');
setfig('FontSize', 12, 'LineWidth', 1);
grid;
printpdf('radiatterer_pf');

figure;
formatfig([14, 4], [0 .5 .5 .5]);
plot(freqvec(1:2000), P_conv, freqvec, P_fmm/4);
legend({'conv BEM', 'FMM'}, 'Location', 'NorthWest');
xlabel('Frequency [Hz]');
ylabel('Radiated power');
setfig('FontSize', 12, 'LineWidth', 1);
grid;
printpdf('radiatterer_power');

figure;
formatfig([14, 4], [0 .5 .5 .5]);
plot(freqvec(1:2000), iters_conv, freqvec, iters_fmm);
legend({'conv BEM', 'FMM'}, 'Location', 'NorthWest');
xlabel('Frequency [Hz]');
ylabel('iterations');
setfig('FontSize', 12, 'LineWidth', 1);
grid;
printpdf('radiatterer_iters');
