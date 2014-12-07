clear;

dim = 3;

source = create_sphere_boundary(1,15);
receiver = translate_mesh(source, [1 1 1]);

y = source.Nodes(:,2:4);
x = receiver.Nodes(:,2:4);

nLeaf = 100;
RowTree = build_cluster_tree(x, nLeaf, 'oc');
ColTree = build_cluster_tree(y, nLeaf, 'oc');

Admit = @(C1, C2)is_admissible_bb(C1, C2, .8);
[B_near, B_far] = build_dual_block_tree(RowTree, ColTree, Admit);

cExp = 5;
rExp = 4;
y0 = bb_tree_cheb_nodes(ColTree, cExp, dim);
x0 = bb_tree_cheb_nodes(RowTree, rExp, dim);

M2L = bb_M2L(B_far, x0, y0, rExp, cExp, @laplace_kernel);
M2M = bb_M2M(ColTree, y0, cExp);
