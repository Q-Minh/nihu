clear;

fdim = 2;

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
    'fvar', eye(fdim), 'svar', 0.0625*[2 1; 1 4], 'cheb_order', 5);
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
n_modes = 100;
fprintf('Computing %d eigenvalues ... ', n_modes); tic;
Afun = @(x)mex_fun('mvp', x);
opts = struct('issym', true);
[phi, lam] = eigs(Afun, size(W,1), W, n_modes, 'lm', opts);
fprintf('Ready in %.3f seconds.\n', toc);

%// Sort the eigenvalues
[lam, i] = sort(diag(lam), 'descend');
phi = phi(:, i);

%% Plot eigenvalues
figure;
plot(lam);
xlabel('Index');
ylabel('Eigenvalue \lambda');
%% Plot a mode
figure;
for n = 1 : 8
    subplot(2, 4, n)

p = reshape(phi(:,n), fdim, []).';
P = sqrt(dot(p, p, 2));
plot_mesh(surface, P);
hold on;
c = centnorm(surface);
quiver(c(:,1), c(:,2), p(:,1), p(:,2), 'Color', 'k');
shading flat;
axis equal off
end
%% Realize
n_real = 8;
xi = phi * diag(sqrt(lam)) * randn(n_modes, n_real);
figure;
c = centnorm(surface);
for n = 1 : n_real
    subplot(2, 4, n)
    p = reshape(xi(:,n), fdim, []).';
P = sqrt(dot(p, p, 2));
plot_mesh(surface, P);
hold on;

quiver(c(:,1), c(:,2), p(:,1), p(:,2), 'Color', 'k');
shading flat;
axis equal off

end
%% Clean up
mex_fun('cleanup');
