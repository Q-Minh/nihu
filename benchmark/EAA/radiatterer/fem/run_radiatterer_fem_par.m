%%
clear;
%% Parameters
% Geometry
L = [2.5 2.0 1.7];          % Size of radiatterer (x,y,z) [m]
Ladd = .4;                  % Additional size for external space
center = L/2;               % Center position
% Mesh
Le = 10e-2 / 2;             % Element size 
tol = Le/100;               % Tolerance for meshing
ie_P = 6;                   % Infinite element order
% Material
rho = 1.3;                  % Density [kg/m3]
c = 340;                    % Speed of sound [m/s]
% Solution parameters
f_min = .5;
f_max = 3;
f_step = .5;
method = 'direct';
n_proc = 4;
%% Create FE model
fprintf('Creating the FE model ... '); tic;
[mesh, points] = create_radiatterer_fem(Le, Ladd);

% Get the outer boundary of the brick
b = mesh_section(get_boundary(mesh), ...
    [[0 0 0]-Ladd + tol; L + Ladd- tol], 'none');
b = drop_mesh_IDs(drop_unused_nodes(b));

% Project IE mesh
ie_mesh = project_ie_mesh(b, ie_P, 'point', center);
ie_mesh.Properties = [1 2 ie_P 1 0 0];

% Assemble model
model = merge_coincident_nodes(join_meshes(mesh, ie_mesh), tol);
model = drop_mesh_IDs(model);

model.Materials = [1 1 rho c 0 0];

% Assemble surface mesh
surf_mesh = mesh_section(get_boundary(model), ...
    [[0 0 0]-Ladd + tol; L + Ladd- tol], 'all');
% Node indices on the surface
[surf_model, surf_idx] = drop_unused_nodes(surf_mesh);
fprintf('Ready. (%.2f s)\n', toc);

save('data/model', 'model', 'surf_model');

%% Assemble field interpolation
fprintf('Assembling field interpolation ...\n'); tic;
[El, Xi] = fe_invmap(model, points, 1e-3);
A_int = fe_interp(model, El, Xi); 
fprintf('Ready (%.2f s)\n', toc);

%% Assemble system matrices
fprintf('Assembling system matrices ...\n'); tic;
[M, K, C, DOF] = model_mk(model);
A = model_a(surf_mesh);
fprintf('Ready (%.2f s)\n', toc);

%% Solve
n_nodes = size(model.Nodes, 1);
v = zeros(n_nodes, 1);
v(surf_idx) = 1e-3;

freq_vec = f_min : f_step : f_max;
n_freq = length(freq_vec);
freq_vec = [freq_vec, nan(1, mod(abs(n_freq - n_proc), n_proc))];

freq_mat = reshape(freq_vec, n_proc, []).';

n_surf = size(surf_idx, 1);
n_field = size(points, 1);

%% Preallocate
proc_freqs = cell(n_proc, 1);
p_surf = cell(n_proc, 1);
p_field = cell(n_proc, 1);
flags = cell(n_proc, 1);
iters = cell(n_proc, 1);
relres = cell(n_proc, 1);
soltimes = cell(n_proc, 1);
%% Open pool
clear_pool();
open_pool(n_proc);
%% Solve
parfor i_proc = 1 : n_proc
    proc_freqs{i_proc} = freq_mat(:, i_proc);
    proc_freqs{i_proc} = proc_freqs{i_proc}(~isnan(proc_freqs{i_proc}));
    proc_n_freqs = length(proc_freqs{i_proc});
    p_surf{i_proc} = nan(n_surf, proc_n_freqs);
    p_field{i_proc} = nan(n_field, proc_n_freqs);
    flags{i_proc} = nan(1, proc_n_freqs);
    iters{i_proc} = nan(1, proc_n_freqs);
    relres{i_proc} = nan(1, proc_n_freqs);
    soltimes{i_proc} = nan(2, proc_n_freqs);
    
    for i_freq = 1 : proc_n_freqs
        if (i_proc == 1)
            fprintf('%4d / %4d, f = %.1f Hz ', i_freq, proc_n_freqs, proc_freqs{i_proc}(i_freq));
        end
        f = proc_freqs{i_proc}(i_freq);
        om = 2*pi*f;
        S = K + 1i*om*C - om^2 * M;
        q = -1i*om*A*v;
        switch method
            case 'direct'
                soltimes{i_proc}(1, i_freq) = 0;
                tic;
                p_sol = S \ q;
                soltimes{i_proc}(2, i_freq) = toc;
            case 'bicgstab'
                % Perform LU
                tic;
                [L,U] = ilu(S,struct('type','nofill','droptol',1e-6));
                soltimes{i_proc}(1, i_freq) = toc;
                % Solve using Bicgstab
                tic;
                [p_sol, flags{i_proc}(i_freq), relres{i_proc}(i_freq), ...
                    iters{i_proc}(i_freq)] = bicgstab(S, q, 1e-6, 3000, L, U);
                soltimes{i_proc}(2, i_freq) = toc;
        end
        % Calculate field pressure
        p_field{i_proc}(:, i_freq) = A_int * p_sol;
        % Calculate and save surface pressure
        p_surf{i_proc}(:, i_freq) = p_sol(surf_idx);
        
        if (i_proc == 1)
            fprintf('Ready. (%.2f s)\n', sum(soltimes{i_proc}(:, i_freq)));
        end
   end
end

%% Finish
clear_pool();

%% Collect data
p_s = zeros(n_surf, n_freq);
p_f = zeros(n_field, n_freq);
t_sol = zeros(2, n_freq);
for i_proc = 1 : n_proc
    [~, idx] = ismember(freq_vec, proc_freqs{i_proc});
    p_f(:, idx) = p_field{i_proc};
    t_sol(:, idx) = soltimes{i_proc};
    
    % Save surface pressures
    for i = 1 : length(proc_freqs{i_proc})
        f = proc_freqs{i_proc}(i);
        p_s = p_surf{i_proc}(:, i);
        save(sprintf('data/p_surf_%05d', round(f*10)), 'p_s', 'f');
    end
end

%%
%figure;
%plot(freqs, p_dB);