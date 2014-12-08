clear;

dim = 3;

display = true;

%%
source = create_sphere_boundary(1,30);
receiver = create_brick_boundary([1 1 1], [40 40 40]);
receiver = translate_mesh(receiver, [.5 .5 .5]);
y = centnorm(source);
x = centnorm(receiver);

%%
nLeaf = 50;
RowTree = build_cluster_tree(x, nLeaf, 'oc');
ColTree = build_cluster_tree(y, nLeaf, 'oc');

Admit = @(C1, C2)is_admissible_bb(C1, C2, .8);
[B_near, B_far] = build_dual_block_tree(RowTree, ColTree, Admit);

cExp = 5;
rExp = 4;
y0 = bb_tree_cheb_nodes(ColTree, cExp, dim);
x0 = bb_tree_cheb_nodes(RowTree, rExp, dim);

kernel = @laplace_kernel;

%% Compute bbFMM sparse matrices
P2P = bb_P2P(B_near, RowTree, ColTree, x, y, kernel);
P2M = bb_P2M(y, ColTree, y0, cExp);
M2M = bb_M2M(ColTree, y0, cExp);
M2L = bb_M2L(B_far, x0, y0, rExp, cExp, kernel);
L2L = bb_M2M(RowTree, x0, rExp);
L2L = L2L.';
L2P = bb_P2M(x, RowTree, x0, rExp);
L2P = L2P.';

%% matrix-vector product
sigma = ones(size(y,1),1);
t0 = tic;
resp = bb_matvec(P2P, P2M, M2L, M2M, L2L, L2P, sigma);
tProd = toc(t0);

%% full product and error
t0 = tic;
resp0 = kernel(x, y) * sigma;
tfull = toc(t0);
error = log10(abs(resp./resp0-1));

%%
if display
    figure;
    plot_mesh(source);
    plot_mesh(receiver, resp);
    light;
    lighting phong;
    set(findall(gcf, 'Type', 'patch'), 'LineStyle', 'none');
end
