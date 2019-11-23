clear;
%close all;

% Le = 2e-2;
% freq = 2000;

Le = 5e-2;
type = 'const';
freq = 600;

fname = sprintf('data/%s_%03dmm/%s_%03dmm_%gHz_pf.res', ...
    type, floor(1000*Le), type, floor(1000*Le), freq);
fid = fopen(fname, 'rt');
header = fscanf(fid, '%g', 2);
data = fscanf(fid, '%g', [2 header(2)]);
pf = complex(data(1,:), data(2,:));
fclose(fid);

fname = sprintf('data/%s_%03dmm/%s_%03dmm_%gHz_ps.res', ...
    type, floor(1000*Le), type, floor(1000*Le), freq);
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
switch type
    case 'const'
        ps_plot = ps;
    case 'gauss'
        ps_plot = mean(reshape(ps, 4, []),1).';
end

fig = figure;
formatfig(fig, [14 9]);
plot_mesh(radiatterer, 20*log10(abs(ps_plot(:)/2e-5))); 
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

%print -dpng radiatterer_5000 -r600

