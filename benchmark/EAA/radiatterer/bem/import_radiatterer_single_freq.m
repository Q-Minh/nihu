clear;

%% import mesh and field
mesh = import_off_mesh('data/radiatterer_1cm_quad.off');
field = import_off_mesh('data/radiatterer_plane_1cm_quad.off');

%% import pressures
ps = import_response('data/radiatterer_single_freq_4000ps.res');
pf = import_response('data/radiatterer_single_freq_4000pf.res');

%% plot result
fig = figure;
formatfig(fig, [13 7], [0 .5 .7 .5]);
hold on;
plot_mesh(mesh, 20*log10(abs(ps/2e-5)));
plot_mesh(field, 20*log10(abs(pf/2e-5)));
axis equal;
shading flat;
cb = colorbar;
ylabel(cb, 'SPL [dB]');
caxis([80 110]);
light
lighting phong
view([55 25]);
axis off;
print -dpng -r600 radiatterer_4k