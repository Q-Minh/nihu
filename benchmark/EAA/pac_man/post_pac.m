clear;

%% Parameters
Le = 3e-3;
freq = 4000;
c = 343.21;
type = 'line';
field = 'field';

% Prefix for output files
pattern = sprintf('pac_man_%03dmm_%s_%05dHz', Le*1000, type, freq);
surf_mesh_name = fullfile('mesh', sprintf('pac_man_surf_%03dmm.off', Le*1000));
result_name = fullfile('data', sprintf('%s_%s_result.mat', pattern, field));

%% Mesh import
% For the line problem, the field mesh is also needed
switch field
    case 'field'
        field_mesh_name = fullfile('mesh', sprintf('pac_man_%s_%03dmm_quad.off', field, Le*1000));
    case 'dir_field'
        field_mesh_name = fullfile('mesh', sprintf('pac_man_%s_%03dmm.off', field, Le*1000));
    otherwise
        error('NiHu');
end
fprintf('Importing surface and field OFF meshes ...'); tic;
surf_mesh = import_off_mesh(surf_mesh_name);
field_mesh = import_off_mesh(field_mesh_name);
fprintf('Ready in %.2f seconds.\n', toc);

%% Result import
load(result_name);

%% Compute incident
switch type
    case 'rad'
        % Incident field is zero
        % Nothing to do
    case 'line'
        % Compute incident field
        k = 2*pi*freq/c;
        [~, pf_inc] = create_pac_man_exc(surf_mesh, field_mesh, k, type);
        pf_inc = pf_inc * (4j);
        pf_scat = pf + pf_inc;
end

%% plot scattered field
fig = figure;
formatfig(fig, [9 7], [0 .5 .5 .5]);
hold on;
plot_mesh(surf_mesh);
plot_mesh(field_mesh, 20*log10(abs(pf/2e-5)));
axis equal;
shading flat;
cb = colorbar;
ylabel(cb, 'SPL [dB]');
cx = caxis;
%caxis(max(cx) + [-40 0]);
caxis([35 75]);
axis off;
set(gca, 'FontSize', 8);
%print('-dpng', '-r600', sprintf('%s_abs.png', pattern));
%close(fig);


%% timing

% fig = figure;
% formatfig(fig);
% bar(levels, times);
% set(gca, 'yScale', 'log');
% legend({'M2M', 'M2L', 'L2L'}, 'Location', 'NorthWest');
% title(sprintf('Sum: %.2g s', sum(sum(times))));
% xlabel('level');
% ylabel('t [s]');
% 
% printpdf(sprintf('%s_rad_times', pattern));



% %% load and plot scattered field
% load(sprintf('%s_line', pattern), 'result_line', 'ps_line', 'pf_line', 'pac', 'field');
% [pac, field, k, q_surf, qs_scat_line, pf_in_line] = create_pac_man(Le);
% 
% fig = figure;
% formatfig(fig, [12 8], [0 .5 .5 .5]);
% hold on;
% plot_mesh(pac);
% plot_mesh(field, 20*log10(abs((pf_line + pf_in_line)/2e-5)));
% axis equal;
% shading flat;
% cb = colorbar;
% ylabel(cb, 'SPL [dB]');
% cx = caxis;
% caxis(max(cx) + [-40 0]);
% axis off;
% print('-dpng', '-r600', sprintf('%s_line_abs', pattern));
% % close(fig);
% 
% %% timing
% txt = textscan(result_line, '%s', 'delimiter', '\n');
% txt = txt{1};
% rows = strncmp('Level #', txt, length('Level #'));
% from = find(diff(rows) == 1, 1, 'last');
% to = find(diff(rows) == -1, 1, 'last');
% txt = txt(from+1 : to+1, :);
% 
% for k = 1 : length(txt) - 1
%     data(k,:) = sscanf(txt{k}, 'Level #%d: M2M (%d), M2L (%d), L2L (%d)');
% end
% 
% levels = data(3:end,1);
% times = data(3:end,2:end);
% times = times / 4 / 1e6;
% 
% fig = figure;
% formatfig(fig);
% bar(levels, times);
% set(gca, 'yScale', 'log');
% legend({'M2M', 'M2L', 'L2L'}, 'Location', 'NorthWest');
% title(sprintf('Sum: %.2g s', sum(sum(times))));
% xlabel('level');
% ylabel('t [s]');
% setfig(fig, 'FontSize', 14);
% grid;
% 
% printpdf(sprintf('%s_line_times', pattern));
