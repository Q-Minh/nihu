clear;

%%
directory = 'data';
meshname = 'radiatterer_10cm_quad.off';
pattern = 'gauss_quad_bm_10cm';

%% read mesh
mesh = import_off_mesh(fullfile(directory, meshname));

%% read field point pressures
files = dir(fullfile(directory, sprintf('%s_*pf.res', pattern)));
for f = 1 : length(files)
    progbar(1, length(files), f);
    fname = files(f).name;
    data = importdata(fullfile(directory, fname));
    frqs_pf(f) = data(1);
    pf(:,f) = complex(data(2:2:end), data(3:2:end));
end

[freqs_pf, i] = sort(frqs_pf);
pf = pf(:,i);

%% read surface pressure
files = dir(fullfile(directory, sprintf('%s_*ps.res', pattern)));
for f = 1 : length(files)
    progbar(1, length(files), f);
    fname = files(f).name;
    data = importdata(fullfile(directory, fname));
    frqs_ps(f) = data(1,1);
    iters(f) = data(1,2);
    ps(:,f) = complex(data(2:end,1), data(2:end,2));
end

[freqs_ps, i] = sort(frqs_ps);
ps = ps(:,i);
iters = iters(i);

%%
save(fullfile('data', sprintf('results_%s', pattern)), 'freqs_pf', 'pf', 'freqs_ps', 'ps');

%%
figure(1);
hold on;
formatfig();
plot(freqs_pf, 20*log10(abs(pf)/2e-5), '-');
setfig('LineWidth', 1);
grid;
% legend(num2str((1:size(pf,1))'));
xlabel('Frequency [Hz]');
ylabel('Sound pressure level [dB]');
legend(num2str((1:size(pf,1))'));

%%
[~, idx] = max(abs(pf(1,:)));
f = freqs_pf(idx);
f = 79;
[~, idx] = min(abs(freqs_ps - f));
figure;
plot_mesh(mesh, 20*log10(abs(ps(:,idx)/2e-5)));
shading flat;
plot_mesh(surface2wireframe(mesh, pi/6));
axis equal tight;
