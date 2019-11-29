clear;

%// create a slab
surface = create_slab([1, 1], [150, 150]);
surface = drop_mesh_IDs(drop_unused_nodes( ...
    mesh_section(surface, [.5 .3 -inf; inf inf inf], 'none')));

%// Extract C++ format mesh with volume elements 
[r_nodes, r_elems] = extract_core_mesh(surface, 'volume');
%% Initialize and set the parameters
mex_fun = @covariance_2d_bbfmm_mex;
mex_fun('init');
mex_fun('set', ...
    'sigma', 1, 'cov_length', 0.25, 'cheb_order', 5);
%% Set the mesh
mex_fun('mesh', r_nodes, r_elems);
%% Create the cluster tree
tree_depth = 7;
fprintf('Creating cluster tree ... '); tic;
mex_fun('tree', 'divide_depth', tree_depth);
fprintf('Ready in %.3f seconds\n', toc);
mex_fun('print_tree');
%% Assemble the FMM matrix with operator precomputation
fprintf('Assembling FMM matrix ... '); tic;
mex_fun('matrix');
fprintf('Ready in %.3f seconds\n', toc);
mex_fun('print_times');
W = mex_fun('get_sparse_identity');
%% Compute eigenvalues using Matlab's eigs
n_modes = 50;
fprintf('Computing %d eigenvalues ... ', n_modes); tic;
Afun = @(x)mex_fun('mvp', x);
[phi, lam] = eigs(Afun, size(W,1), W, n_modes, 'lm');
fprintf('Ready in %.3f seconds.\n', toc);

%// Sort the eigenvalues
[lam, i] = sort(diag(lam), 'descend');
phi = phi(:, i);

%% Plot
figure;
plot(lam);
xlabel('Index');
ylabel('Eigenvalue \lambda');

figure;
plot_mesh(surface, phi(:, 24));
shading flat;
%% Clean up
mex_fun('cleanup');