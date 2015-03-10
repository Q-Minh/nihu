clear;

dim = 3;

display = true;

%%
source = create_sphere_boundary(1,30);
x = centnorm(source);

%%
nLeaf = 35;
CTree = build_cluster_tree(x, nLeaf, 'oc');
Admit = @(C1, C2)is_admissible_bb(C1, C2, 1.1);
[B_near, B_far] = build_dual_block_tree(CTree, CTree, Admit);
nExp = 3;
x0 = bb_tree_cheb_nodes(CTree, nExp, dim);

%% Compute bbFMM sparse matrices
kernel = @laplace_kernel;

P2P = bb_P2P(B_near, CTree, CTree, x, x, kernel);
P2M = bb_P2M(x, CTree, x0, nExp);
M2M = bb_M2M(CTree, x0, nExp);
M2L = bb_M2L(B_far, x0, x0, nExp, nExp, kernel);

%% matrix-vector product
sigma = ones(size(x,1),1);
t0 = tic;
resp = bb_matvec(P2P, P2M, M2L, M2M.', L2L, P2M.', sigma);
tProd = toc(t0);

%% Solve back
Afun = @(sigma)bb_matvec(P2P, P2M, M2L, M2M, L2L, L2P, sigma);
t0 = tic;
[sigmaback, flag, relres, iter, resvec] = gmres(Afun, resp);
tSol = toc(t0);

% %% full product and error
% t0 = tic;
% resp0 = kernel(x, x) * sigma;
% tFull = toc(t0);
% error = log10(abs(resp./resp0-1));

%%
if display
    figure;
    plot_mesh(source, resp);
    light;
    lighting phong;
    set(findall(gcf, 'Type', 'patch'), 'LineStyle', 'none');
end
