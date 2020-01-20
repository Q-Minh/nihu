% run_pac.m
% Run the Pac-Man at a selected frequency and excitation case
% Requires that mesh is generated 
clear;

Le = 3e-3;
freq = 2000;
field = 'dir_field';
type = 'line';

if isunix
    runner_name = 'sh runner_linux.sh';
    n_threads = 32;
else
    runner_name = 'runner_msvc.bat';
end


% Prefix for output files
pattern = sprintf('pac_man_%03dmm_%s_%05dHz', Le*1000, type, freq);
surf_mesh_name = fullfile('mesh', sprintf('pac_man_surf_%03dmm.off', Le*1000));
surf_mesh = import_off_mesh(surf_mesh_name);
field_mesh_name = fullfile('mesh', sprintf('pac_man_%s_%03dmm.off', field, Le*1000));

%% create excitation
c = 343.21;         % Reference speed of sound for Pac-Man [m/s]

k = 2*pi*freq/c;       % Wave number [rad/m]

q_surf = create_pac_man_exc(surf_mesh, [], k, type);

surf_exc_name = fullfile('data', sprintf('%s.xct', pattern));
export_excitation(q_surf, k, surf_exc_name);

surf_result_name = fullfile('data', sprintf('%s_ps.res', pattern));
field_result_name = fullfile('data', sprintf('%s_%s_pf.res', pattern, field));

%% call the runner that solves the problem
stdout_file_name = fullfile('data', sprintf('%s_%s_stdout.txt', pattern, field));
if (isunix)
command = sprintf('%s %s %s %s %s %s %d > %s', runner_name, ...
    surf_mesh_name, surf_exc_name, surf_result_name, ...
    field_mesh_name, field_result_name, n_threads, stdout_file_name);
else
    command = sprintf('%s %s %s %s %s %s > %s', runner_name, ...
        surf_mesh_name, surf_exc_name, surf_result_name, ...
        field_mesh_name, field_result_name, stdout_file_name);
end
fprintf('Calling 2D WB FMM executable ...\n'); tic;
[status, result] = system(command);
fprintf('Completed in %.2f seconds.\n', toc);
if status == 0 && isempty(result)
    fprintf('Command succeeded.\n');
else
    fprintf('Error occured during execution, result: %s\n', result);
end

%% import the results to convert to Matlab
ps = import_response(surf_result_name);
pf = import_response(field_result_name);

%% Process timing info from standard output
fid = fopen(stdout_file_name, 'rt');
txt = textscan(fid, '%s', 'delimiter', '\n');
txt = txt{1};
fclose(fid);

% Find number of threads
rows = find(strncmp('Expanding to ', txt, length('Expanding to ')));
if (~isempty(rows))
    n_threads = sscanf(txt{rows(end)}, 'Expanding to %d threads');
else
    n_threads = 1;
end

% Find MVP wallclock and CPU times
mvp_wc_str = 'MVP wall clock time:';
idx = find(strncmp(mvp_wc_str, txt, length(mvp_wc_str)), 1, 'first');
if (~isempty(idx))
    mvp_wc_t = sscanf(txt{idx}, [mvp_wc_str ' %f']);
else
    mvp_wc_t = nan;
end

mvp_cpu_str = 'MVP CPU time:';
idx = find(strncmp(mvp_cpu_str, txt, length(mvp_cpu_str)), 1, 'first');
if (~isempty(idx))
    mvp_cpu_t = sscanf(txt{idx}, [mvp_cpu_str ' %f']);
else
    mvp_cpu_t = nan;
end

n_src_str = '#Source Nodes:';
idx = find(strncmpi(n_src_str, txt, length(n_src_str)), 1, 'last');
if (~isempty(idx))
    n_src = sscanf(txt{idx}, [n_src_str ' %d']);
else
    n_src = nan;
end

n_rec_str = '#Receiver Nodes:';
idx = find(strncmpi(n_src_str, txt, length(n_src_str)), 1, 'last');
if (~isempty(idx))
    n_rec = sscanf(txt{idx}, [n_src_str ' %d']);
else
    n_rec = nan;
end

% Find timing info for levels
rows = strncmp('Level #', txt, length('Level #'));
from = find(diff(rows) == 1, 1, 'last');
to = find(diff(rows) == -1, 1, 'last');
txt = txt(from+1 : to+1, :);

data = zeros(length(txt)-1, 4);
for k = 1 : length(txt) - 1
    data(k,:) = sscanf(txt{k}, 'Level #%d: M2M (%d), M2L (%d), L2L (%d)');
end

levels = data(3:end,1);
times = data(3:end,2:end);
% Note: should number of threads be taken into account?
% Convert from microseconds to seconds
times = times / 1e6;

%% Save the result
save(fullfile('data', sprintf('%s_%s_result', pattern, field)), ...
    'ps', 'pf', 'levels', 'times', 'n_threads', 'mvp_cpu_t', 'mvp_wc_t');

%% Save timing data separately
td_name = fullfile('timing_data', sprintf('%s_%02d_timing', pattern, n_threads));
save(td_name, 'levels', 'times', 'n_threads', 'mvp_cpu_t', 'mvp_wc_t', 'n_src', 'n_rec');
