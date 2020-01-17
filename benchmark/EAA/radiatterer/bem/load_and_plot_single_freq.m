clear;

Le_s = 6.25e-3;         % Element size on surface [m]
Le_f = 3.125e-3;        % Element size on field [m]
type = 'const';         % Const / Gauss elements
freq = 6000;            % Frequency [Hz]
method = 'fmm';         % FMM / Conv. method

% Import surface mesh
fprintf('Importing surface mesh ... '); tic;
meshname = sprintf('mesh/radiatterer_%03dmm_quad.off', floor(1000*Le_s));
radiatterer = import_off_mesh(meshname);
wf = surface2wireframe(radiatterer, pi/3);
fprintf('Ready in %.2f seconds.\n', toc);

% Read surface result
fname = sprintf('data_%s/%s_%03dmm/%s_%03dmm_%gHz_ps.res', ...
    method, type, floor(1000*Le_s), type, floor(1000*Le_s), freq);
fprintf('Loading surface result ... '); tic;
fid = fopen(fname, 'rt');
header = fscanf(fid, '%g', 2);
data = fscanf(fid, '%g', [2 header(2)]);
ps = complex(data(1,:), data(2,:));
fclose(fid);
fprintf('Ready in %.2f seconds.\n', toc);

% Load the field for FMM
if (strcmp(method, 'fmm'))
    % Import field mesh
    fprintf('Importing field mesh ... '); tic;
    fieldname = sprintf('mesh/radi_plane_%03dmm_quad.off', floor(1000*Le_f));
    field = import_off_mesh(fieldname);
    fprintf('Ready in %.2f seconds.\n', toc);
    
    % Read field result
    fname = sprintf('data_%s/%s_%03dmm/%s_%03dmm_%gHz_pf.res', ...
        method, type, floor(1000*Le_s), type, floor(1000*Le_f), freq);
    fprintf('Loading field result ... '); tic;
    fid = fopen(fname, 'rt');
    header = fscanf(fid, '%g', 2);
    data = fscanf(fid, '%g', [2 header(2)]);
    pf = complex(data(1,:), data(2,:));
    fclose(fid);
    fprintf('Ready in %.2f seconds.\n', toc);
end

%% Rearrange surface pressure for plotting
switch type
    case 'const'
        ps_plot = ps;
    case 'gauss'
        ps_plot = mean(reshape(ps, 4, []),1).';
end

%%
fig = figure;
formatfig(fig, [14 9]);
plot_mesh(radiatterer, 20*log10(abs(ps_plot(:)/2e-5))); 
plot_mesh(wf); 
if (strcmp(method, 'fmm'))
    plot_mesh(field, 20*log10(abs(pf(:)/2e-5))); 
end
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

