clear;

Le = 2e-3;
freq = 250;
c = 343.21;
type = 'line';
field = 'dir_field';

% Prefix for output files
pattern = sprintf('pac_man_%03dmm_%s_%05dHz', Le*1000, type, freq);
surf_mesh_name = fullfile('mesh', sprintf('pac_man_surf_%03dmm.off', Le*1000));
result_name = fullfile('data', sprintf('%s_%s_result.mat', pattern, field));

%% Mesh import
% For the line problem, the field mesh is also needed
switch field
    case 'field'
        field_mesh_name = fullfile('mesh', sprintf('pac_man_%s_%03dmm_quad.off', field, Le*1000));
    case {'dir_field', 'ref_field'}
        field_mesh_name = fullfile('mesh', sprintf('pac_man_%s_%03dmm.off', field, Le*1000));
    otherwise
        error('NiHu');
end
fprintf('Importing surface and field OFF meshes ... '); tic;
surf_mesh = import_off_mesh(surf_mesh_name);
field_mesh = import_off_mesh(field_mesh_name);
fprintf('Ready in %.2f seconds\n', toc);

fprintf('Importing field point result ... '); tic;
load(result_name, 'pf');
fprintf('Ready in %.2f seconds\n', toc);

switch type
    case 'rad'
        % Incident field is zero
        % Nothing to do
    case 'line'
        % Compute incident field
        k = 2*pi*freq/c;
        [~, pf_inc] = create_pac_man_exc(surf_mesh, field_mesh, k, type);
        pf_inc = pf_inc * (4j);
        %pf = pf + pf_inc;
end

%% Load database
switch type
    case 'rad'
        db_name = fullfile('ref_impl', 'reference_solution_surface_vibration.csv');
    case 'line'
        x_ref = load(fullfile('ref_impl', 'python', 'data', 'nodes.txt'), '-ascii');
        p_ref_name = fullfile('ref_impl', 'python', 'data', sprintf('freq_%.1f.txt', freq));
        fid = fopen(p_ref_name, 'rt');
        txt = textscan(fid, '%s');
        txt = txt{1};
        fclose(fid);
        txt = cellfun(@strtrim, txt, 'UniformOutput', false);
        txt = cellfun(@(x)x(2:end-1), txt, 'UniformOutput', false);
        p_ref = cellfun(@str2double,txt);
    otherwise
end

phi_ref = atan2(x_ref(:,2) ,x_ref(:,1));

%% Plot comparison
field_xyz = centnorm(field_mesh);
phi_field = atan2(field_xyz(:,2), field_xyz(:,1));
[~, idx] = sort(phi_field);

figure;
plot(phi_field(idx), real(pf(idx)), '-');
hold on
plot(phi_ref, real(p_ref), '.');
legend({'NiHu solution', 'Reference solution'}, 'Location', 'SouthWest');

%% Plot polar
figure;
formatfig([5 5], [.5 .5 .5 .5]);
polarplot(phi_field(idx), abs(pf(idx))/max(abs(pf(idx))),...
    'LineWidth', .75);
hold on 
step = 15;
polarplot(phi_ref(1:step:end), abs(p_ref(1:step:end))/max(abs(p_ref)), 'o', ...
    'MarkerSize', 2, 'Color', 'black');
set(gca, 'thetatick', 0 : 45 : 360);
set(gca, 'rtick', 0:.2:1);
set(gca, 'rticklabel', []);
text(-pi/3, 1.2, sprintf('%d Hz', freq), 'FontSize', 8);
set(gca, 'FontSize', 8);
printpdf(sprintf('pac_compare_%04d_Hz', freq));
