%TEST_DUAL_3D test ACA with different source and receiver 3D meshes

clear;

%%
n = 50;
source = create_circle_boundary(1, n);
source = extrude_mesh(source, [0 0 2*pi/n], n);
receiver = scale_mesh(source, [.5 .5 1]);
xs = centnorm(source);
xr = centnorm(receiver);
Ns = size(xs,1);
Nr = size(xr,1);

%% build cluster trees
nLeaf = 20;
t0 = tic;
SourceTree = build_cluster_tree(xs, nLeaf);
[SourceClusters, SourcePerm] = sort_cluster_tree(SourceTree);
ReceiverTree = build_cluster_tree(xr, nLeaf);
[ReceiverClusters, ReceiverPerm] = sort_cluster_tree(ReceiverTree);
tTree = toc(t0);
fprintf('%.3g s needed to build cluster tree\n', tTree);

%% build block tree
t0 = tic;
[B_near, B_far] = build_dual_block_tree(ReceiverTree, SourceTree);
tBtree = toc(t0);
fprintf('%.3g s needed to build block tree\n', tBtree);

%% visualise far field block tree
display_block_structure(ReceiverClusters, SourceClusters, B_far);
set(gcf, 'renderer', 'painters');

%%
k = min(mesh_kmax(source));
M = @(i,j)helmholtz_matrix(xr(i,:),xs(j,:),k);

exc = ones(Ns,1);
resp = zeros(Nr, 1);

t0 = tic;
for b = 1 : size(B_near,1)
    I = ReceiverTree(B_near(b,1)).ind;
    J = SourceTree(B_near(b,2)).ind;
    resp(I) = resp(I) + M(I,J) * exc(J);
end
tNear = toc(t0);
fprintf('%.3g s needed to compute near field\n', tNear);

eps = 1e-2;
t0 = tic;
for b = 1 : size(B_far,1)
    progbar(1, size(B_far,1), b);
    I = ReceiverTree(B_far(b,1)).ind;
    J = SourceTree(B_far(b,2)).ind;
    [U, V] = lowrank_approx_block(M, I, J, eps);
    resp(I) = resp(I) + U * (V' * exc(J));
end
tAca = toc(t0);
fprintf('%.3g s needed to compute far field\n', tAca);

%%
t0 = tic;
resp0 = M(1:Nr,1:Ns) * exc;
tFull = toc(t0);
fprintf('%.3g s needed to compute full product\n', tFull);

logeps = log10(abs(resp./resp0-1));

figure;
plot_mesh(receiver, logeps);
c = colorbar;
ylabel(c, 'log10 error [-]');
caxis([-8 0]);
set(gcf, 'renderer', 'painters');
