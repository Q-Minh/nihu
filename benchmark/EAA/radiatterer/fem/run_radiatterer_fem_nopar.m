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
ie_P = 4;                   % Infinite element order
% Material
rho = 1.3;                  % Density [kg/m3]
c = 340;                    % Speed of sound [m/s]
% Solution parameters
f_min = .5;
f_max = 1000;
f_step = .5;
method = 'direct';

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

freqs = f_min : f_step : f_max;
n_freq = length(freqs);

n_surf = size(surf_idx, 1);
n_field = size(points, 1);

%% Preallocate
p_field = nan(n_field, n_freq);
soltimes = nan(2, n_freq);
iters = nan(1, n_freq);
relres = nan(1, n_freq);
flags = nan(1, n_freq);

%% Solve
for i_freq = 1 : n_freq
    f = freqs(i_freq);
    fprintf('%4d / %4d, f = %.1f Hz ', i_freq, n_freq, f);
   
    om = 2*pi*f;
    S = K + 1i*om*C - om^2 * M;
    q = -1i*om*A*v;
    switch method
        case 'direct'
            soltimes(1, i_freq) = 0;
            tic;
            p_sol = S \ q;
            soltimes(2, i_freq) = toc;
        case 'bicgstab'
            % Perform LU
            tic;
            [L,U] = ilu(S,struct('type','nofill','droptol',1e-6));
            soltimes(1, i_freq) = toc;
            % Solve using Bicgstab
            tic;
            [p_sol, flags(i_freq), relres(i_freq), ...
                iters(i_freq)] = bicgstab(S, q, 1e-6, 10000, L, U);
            soltimes(2, i_freq) = toc;
        case 'gmres'
            % Perform LU
            tic;
            [L,U] = ilu(S,struct('type','nofill','droptol',1e-6));
            soltimes(1, i_freq) = toc;
            % Solve using Bicgstab
            tic;
            inner = 500;
            maxit = 20;
            [p_sol, flags(i_freq), relres(i_freq), ...
                it] = gmres(S, q, inner, 1e-6, maxit, L, U);
            iters(i_freq) = (it(1) - 1)*inner + it(2);
            soltimes(2, i_freq) = toc;
    end
    % Calculate field pressure
    p_field(:, i_freq) = A_int * p_sol;
    % Calculate and save surface pressure
    p_surf = p_sol(surf_idx);
    save(sprintf('data/p_surf_%05d', round(f*10)), 'p_surf');
    fprintf('Ready. (%.2f s)\n', sum(soltimes(:,i_freq)));
end

%% Save the results
save('data/result', 'points', 'freqs', 'p_field', ...
    'method', 'relres', 'iters', 'flags', 'soltimes');

%% Create figure
figure;
p_dB = 20*log10(abs(p_field)/2e-5);
plot(freqs, p_dB, 'LineWidth', 1.5);
xlabel('Frequency [Hz]');
ylabel('p [dBSPL]');
legend(strcat({'point '}, num2str((1:n_field).')));

%% Load the results
clear;
load('data/result');
load('data/model', 'surf_model');

n_surf = size(surf_model.Elements, 1);
[xg, ~, wg, El] = geo2gauss(surf_model, 2);
xi = gaussquad2(2, 4);
Xi = repmat(xi, n_surf, 1);
A = fe_interp(surf_model, El, Xi);

n_freq = length(freqs);
P_rad = zeros(1, n_freq);
for i_freq = 1 : n_freq
    progbar(1, n_freq, i_freq);
    f = freqs(i_freq);
    load(sprintf('data/p_surf_%05d', round(f*10)), 'p_surf');
    
    p_g = A * p_surf;
    % Constant velocity over the surface
    % Negative sign because of normal direction points inwards
    v_const = -1e-3;        
    P_rad(i_freq) = v_const * wg.' * p_g;
end

%% Store 
save('data/result', '-append', 'P_rad');