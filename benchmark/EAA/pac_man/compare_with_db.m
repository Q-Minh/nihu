clear;

Le = 2e-3;
freq = 4000;
c = 343.21;
type = 'rad';
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
        pf = pf + pf_inc;
end

%% Load database
switch type
    case 'rad'
        db_name = fullfile('ref_impl', 'reference_solution_surface_vibration.csv');
    case 'line'
        db_name = fullfile('ref_impl', 'reference_solution_line_source.csv');
    otherwise
end
db_freqs = [NaN NaN 16 31.5 63 125 250 500 1000 2000 4000];
db = csvread(db_name, 1, 0);

phi_db = db(:, 2) * pi / 180;       % Phi as radian
% Correction
phi_db = phi_db + 2.5 * pi /180;
% Shift to -pi ... pi
phi_db = atan2(sin(phi_db), cos(phi_db));
col_db = find(db_freqs == freq, 1, 'first');
p_db = db(:, col_db);
%% Plot comparison
field_xyz = centnorm(field_mesh);
phi_field = atan2(field_xyz(:,2), field_xyz(:,1));
[~, idx] = sort(phi_field);

figure;
plot(phi_field(idx), real(pf(idx)), '-');
hold on
plot(phi_db, real(p_db), '.');
legend({'NiHu solution', 'Reference solution'}, 'Location', 'SouthWest');

%% Plot polar
figure;
polarplot(phi_field(idx), abs(pf(idx)));