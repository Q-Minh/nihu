clear;
close all;
nThrd = 1;
fname = sprintf('timing_data/pac_man_003mm_rad_04000Hz_%02d_timing', nThrd);
load(fname);

fig = figure;
formatfig(fig, [9 5], [1 1 .5 .5]);
bar(levels, times);
set(gca, 'yScale', 'log');
xlabel('level');
ylabel('CPU time [s]');
legend('M2M', 'M2L', 'L2L', 'location', 'northwest');
set(gca, 'FontSize', 8);
grid;
printpdf('pac_level_timing');