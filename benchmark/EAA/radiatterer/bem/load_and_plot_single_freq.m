clear;
close all;

% Le = 2e-2;
% freq = 2000;

Le = .1/15;
freq = 5000;

fname = sprintf('data/const_quad_%03dmm/const_quad_%03dmm_%gHz_pf.res', floor(1000*Le), floor(1000*Le), freq);
fid = fopen(fname, 'rt');
header = fscanf(fid, '%g', 2);
data = fscanf(fid, '%g', [2 header(2)]);
pf = complex(data(1,:), data(2,:));
fclose(fid);

fname = sprintf('data/const_quad_%03dmm/const_quad_%03dmm_%gHz_ps.res', floor(1000*Le), floor(1000*Le), freq);
fid = fopen(fname, 'rt');
header = fscanf(fid, '%g', 2);
data = fscanf(fid, '%g', [2 header(2)]);
ps = complex(data(1,:), data(2,:));
fclose(fid);

meshname = sprintf('data/radiatterer_%03dmm_quad.off', floor(1000*Le));
radiatterer = import_off_mesh(meshname);
wf = surface2wireframe(radiatterer, pi/3);

fieldname = sprintf('data/radi_plane_%03dmm_quad.off', floor(1000*Le));
field = import_off_mesh(fieldname);

%%
fig = figure;
formatfig(fig, [14 9]);
plot_mesh(radiatterer, 20*log10(abs(ps(:)/2e-5))); 
plot_mesh(wf); 
plot_mesh(field, 20*log10(abs(pf(:)/2e-5))); 
view(150,25);
shading flat;
axis equal;
hcb = colorbar;
ylabel(hcb, 'Sound pressure level [dBSPL]');
caxis([-40 0] + 110);
hl = light;
lighting phong;
hl.Position = [2 2 1];

print -dpng radiatterer_5000 -r600

