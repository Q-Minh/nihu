clear;

kernel = @(x,y)laplace_kernel(x,y);

%%
interval = [0 20]';
D = diff(interval); % geometry size
% Nx points in first third
Nx = 1e4;
x = linspace(interval(1,1), interval(2,1)/3, Nx+2)';
x = x(2:end-1); % drop end points
% Ny points in 3rd third
Ny = 1e4;
y = linspace(interval(2)*2/3, interval(2), Ny+2)';
y = y(2:end-1);

nCheb = 5;          % expansion length of Cheby interpolation

%% build cluster tree
nLeaf = 10;         % maximal number of particles on leaf level
depth = round(log2(max(Nx, Ny)/nLeaf));
[tree, m2m] = build_tree(D, depth, nCheb, kernel);

%% build leaf contribution sparse matrices
fprintf(1, 'Computing P2M and L2P...');
[P2M, L2P] = leaf_contribution(x, y, tree(end), nCheb);
fprintf(1, 'done\n');

%%
sigma = rand(Ny,1);
tic;
Multipole = P2M * sigma;
for level = depth : -1 : 2
    tree(level).Local = tree(level).M2L * Multipole;
    if level > 2
        Multipole = reshape(m2m * reshape(Multipole, 2*nCheb, []), [], 1);
    end
end
for level = 3 : depth
    tree(level).Local = tree(level).Local +...
        reshape(m2m' * reshape(tree(level-1).Local, nCheb, []), [], 1);
end
res = L2P * tree(end).Local;
tfmm = toc;

%% Direct method
K = kernel(x, y);
tic;
res0 = kernel(x, y) * sigma;
tanal = toc;

eps = norm(res-res0)/norm(res0);

fprintf(1, 'Error: %g\nTime gain: %g\n', eps, tanal/tfmm);
