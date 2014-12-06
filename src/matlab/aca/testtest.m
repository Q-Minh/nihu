clear;

%%
% if isunix
%     root = '/D';
% else
%     root = 'D:';
% end
% m = import_mesh(fullfile(root, 'research', 'pub',...
%     '2013', 'Boonen2013', 'work', 'industrial', 'data',...
%     'horse.off'));
%%
m = create_sphere_boundary(1, 15);
x = centnorm(m);

%%
t0 = tic;
nLeaf = 25;
Ctree = build_cluster_tree(x, nLeaf);
tCtree = toc(t0);
fprintf('%.3g s needed to build cluster tree\n', tCtree);

%% build block tree
t0 = tic;
[B_near, B_far] = build_block_tree(Ctree);
tBtree = toc(t0);
fprintf('%.3g s needed to build block tree\n', tBtree);

%%
k = min(mesh_kmax(m));
M = @(i,j)helmholtz_matrix(i,j,x,k);

N = size(x,1);
exc = ones(N,1);
resp = zeros(N,1);

t0 = tic;
for i = 1 : size(B_near,1)
    I = Ctree(B_near(i,1)).ind;
    J = Ctree(B_near(i,2)).ind;
    resp(I) = resp(I) + M(I,J) * exc(J);
end
tNear = toc(t0);
fprintf('%.3g s needed to compute near field\n', tNear);

eps = 1e-2;
t0 = tic;
for b = 1 : size(B_far,1)
    progbar(1, size(B_far,1), b);
    I = Ctree(B_far(b,1)).ind;
    J = Ctree(B_far(b,2)).ind;
    [U, V] = lowrank_approx_block(M, I, J, eps);
    resp(I) = resp(I) + U * (V' * exc(J));
end
tACA = toc(t0);
fprintf('%.3g s needed to compute far field\n', tACA);

%%
t0 = tic;
resp0 = M(1:N,1:N) * exc;
tFull = toc(t0);
fprintf('%.3g s needed to compute full field\n', tFull);

logeps = log10(abs(resp./resp0-1));

figure;
plot_mesh(m, logeps);
colorbar;
