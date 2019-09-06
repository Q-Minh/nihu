clear;

%// create a slab
surface = create_slab([1, 1], [300, 300]);
surface = drop_mesh_IDs(drop_unused_nodes( ...
    mesh_section(surface, [.5 .3 -inf; inf inf inf], 'none')));

%// Extract C++ format mesh with volume elements 
[r_nodes, r_elems] = extract_core_mesh(surface, 'volume');
%% Initialize and set the parameters
covariance_2d_bbfmm_matlab('init');
covariance_2d_bbfmm_matlab('set', ...
    'sigma', 1, 'cov_length', 0.25, 'cheb_order', 5);
%% Set the mesh
covariance_2d_bbfmm_matlab('mesh', r_nodes, r_elems);
%% Create the cluster tree
tree_depth = 7;
fprintf('Creating cluster tree ... '); tic;
covariance_2d_bbfmm_matlab('tree', 'divide_depth', tree_depth);
fprintf('Ready in %.3f seconds\n', toc);
%% Assemble the FMM matrix with operator precomputation
fprintf('Assembling FMM matrix ... '); tic;
covariance_2d_bbfmm_matlab('matrix');
fprintf('Ready in %.3f seconds\n', toc);

%% Compute eigenvalues using Matlab's eigs

%// Create diagonal matrix of element sizes
[~, ~, w] = geo2gauss(surface, 1);
W = spdiags(w, 0, size(w,1), size(w,1));

n_modes = 100;
fprintf('Computing %d eigenvalues ... ', n_modes); tic;
Afun = @(x)covariance_2d_bbfmm_matlab('mvp', x);
[phi, lam] = eigs(Afun, size(w,1), W, n_modes, 'lm');
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
covariance_2d_bbfmm_matlab('cleanup');
