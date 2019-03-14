clear;

%%
Le = 3e-3;
pattern = sprintf('data/pac_man_%gmm', 1000*Le);

%% load surface radiation data
load(sprintf('%s_rad_surf', pattern), 'result_rad', 'ps_rad', 'pf_rad', 'pac', 'field');

fig = figure;
formatfig(fig, [12 8], [0 .5 .5 .5]);
hold on;
plot_mesh(pac);
plot_mesh(field, 20*log10(abs(pf_rad/2e-5)));
axis equal;
shading flat;
cb = colorbar;
ylabel(cb, 'SPL [dB]');
cx = caxis;
caxis(max(cx) + [-40 0]);
axis off;
print('-dpng', '-r600', sprintf('%s_rad_abs', pattern));
close(fig);


%% timing
txt = textscan(result_rad, '%s', 'delimiter', '\n');
txt = txt{1};
rows = strncmp('Level #', txt, length('Level #'));
from = find(diff(rows) == 1, 1, 'last');
to = find(diff(rows) == -1, 1, 'last');
txt = txt(from+1 : to+1, :);

for k = 1 : length(txt) - 1
    data(k,:) = sscanf(txt{k}, 'Level #%d: M2M (%d), M2L (%d), L2L (%d)');
end

levels = data(3:end,1);
times = data(3:end,2:end);
times = times / 4 / 1e6;

fig = figure;
formatfig(fig);
bar(levels, times);
set(gca, 'yScale', 'log');
legend({'M2M', 'M2L', 'L2L'}, 'Location', 'NorthWest');
title(sprintf('Sum: %.2g s', sum(sum(times))));
xlabel('level');
ylabel('t [s]');

printpdf(sprintf('%s_rad_times', pattern));

%% load and plot scattered field
load(sprintf('%s_line', pattern), 'result_line', 'ps_line', 'pf_line', 'pac', 'field');
[pac, field, k, q_surf, qs_scat_line, pf_in_line] = create_pac_man(Le);

fig = figure;
formatfig(fig, [12 8], [0 .5 .5 .5]);
hold on;
plot_mesh(pac);
plot_mesh(field, 20*log10(abs((pf_line + pf_in_line)/2e-5)));
axis equal;
shading flat;
cb = colorbar;
ylabel(cb, 'SPL [dB]');
cx = caxis;
caxis(max(cx) + [-40 0]);
axis off;
print('-dpng', '-r600', sprintf('%s_line_abs', pattern));
% close(fig);

%% timing
txt = textscan(result_line, '%s', 'delimiter', '\n');
txt = txt{1};
rows = strncmp('Level #', txt, length('Level #'));
from = find(diff(rows) == 1, 1, 'last');
to = find(diff(rows) == -1, 1, 'last');
txt = txt(from+1 : to+1, :);

for k = 1 : length(txt) - 1
    data(k,:) = sscanf(txt{k}, 'Level #%d: M2M (%d), M2L (%d), L2L (%d)');
end

levels = data(3:end,1);
times = data(3:end,2:end);
times = times / 4 / 1e6;

fig = figure;
formatfig(fig);
bar(levels, times);
set(gca, 'yScale', 'log');
legend({'M2M', 'M2L', 'L2L'}, 'Location', 'NorthWest');
title(sprintf('Sum: %.2g s', sum(sum(times))));
xlabel('level');
ylabel('t [s]');
setfig(fig, 'FontSize', 14);
grid;

printpdf(sprintf('%s_line_times', pattern));
