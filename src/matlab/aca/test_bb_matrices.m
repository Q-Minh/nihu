clear;

% problem dimensionality
dim = 3;

% source and receiver geometries
source = create_sphere_boundary(1,20);
receiver = translate_mesh(source, [2 2 2]);

% extract source and receiver points
y = source.Nodes(:,2:4);
x = receiver.Nodes(:,2:4);

% build oc cluster trees
nLeaf = 100;
RowTree = build_cluster_tree(x, nLeaf, 'oc');
ColTree = build_cluster_tree(y, nLeaf, 'oc');

% build block tree
Admit = @(C1, C2)is_admissible_bb(C1, C2, .8);
[B_near, B_far] = build_dual_block_tree(RowTree, ColTree, Admit);

% get chebyshev nodes for each cluster
cExp = 5;
rExp = 4;
y0 = bb_tree_cheb_nodes(ColTree, cExp, dim);
x0 = bb_tree_cheb_nodes(RowTree, rExp, dim);

% compute bbFMM sparse matrices
P2P = bb_P2P(B_near, RowTree, ColTree, x, y, @laplace_kernel);
M2L = bb_M2L(B_far, x0, y0, rExp, cExp, @laplace_kernel);
M2M = bb_M2M(ColTree, y0, cExp);
P2M = bb_P2M(y, ColTree, y0, cExp);
L2L = bb_M2M(RowTree, x0, rExp);
L2P = bb_P2M(x, RowTree, x0, rExp);
L2L = L2L.';
L2P = L2P.';
