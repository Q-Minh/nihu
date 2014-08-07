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
source = create_line(1, 3000);
receiver = translate_mesh(create_line(1, 3000), [.5 0 0]);
xs = centnorm(source);
xr = centnorm(receiver);

%% build cluster trees
t0 = tic;
SourceTree = build_cluster_tree(xs, 25);
[SourceClusters, SourcePerm] = sort_cluster_tree(SourceTree);
ReceiverTree = build_cluster_tree(xr, 25);
[ReceiverClusters, ReceiverPerm] = sort_cluster_tree(ReceiverTree);
tTree = toc(t0);
fprintf('%.3g s needed to build cluster tree\n', tTree);

%% build block tree
t0 = tic;
[B_near, B_far] = build_block_tree_2(SourceTree, ReceiverTree);
tBtree = toc(t0);
fprintf('%.3g s needed to build block tree\n', tBtree);

%% visualise far field block tree
for b = 1 : size(B_far,1)
    i = ReceiverClusters(B_far(b,1),:);
    j = SourceClusters(B_far(b,2),:);
    patch(j(1) + j(2)*[0 0 1 1]', i(1) + i(2)*[0 1 1 0]', b*ones(4,1));
end

%% visualise near field block tree
v = nan(size(xr,1),1);
w = nan(size(xs,1),1);

sigma = unique(B_near(:,1));
sigma = sigma(end);
v(ReceiverTree(sigma).ind) = 1;
tau = B_near(B_near(:,1) == sigma,2);
for t = 1: length(tau);
    w(SourceTree(tau(t)).ind) = 2;
end

figure;
plot_mesh(source, w);
plot_mesh(receiver, v*100);

shading flat;

%% number of near and far entries
nnear = 0;
for b = 1 : size(B_near,1)
    nnear = nnear + length(ReceiverTree(B_near(b,1)).ind) * length(SourceTree(B_near(b,2)).ind);
end

nfar = 0;
for b = 1 : size(B_far,1)
    nfar = nfar + length(ReceiverTree(B_far(b,1)).ind) * length(SourceTree(B_far(b,2)).ind);
end

%%
k = min(mesh_kmax(m));
M = @(i,j)helmholtz_matrix(i,j,x,k);

N = size(x,1);
exc = ones(N,1);

resp = zeros(N,1);

eps = 1e-2;
t0 = tic;
for tau = 1 : size(B_far,1)
    progbar(1, size(B_far,1), tau);
    I = SourceTree(B_far(tau,1)).ind;
    J = SourceTree(B_far(tau,2)).ind;
    [U, V] = lowrank_approx_block(M, I, J, eps);
    resp(I) = resp(I) + U * (V' * exc(J));
end
tACA = toc(t0);
fprintf('%.3g s needed to compute far field\n', tACA);

%%
t0 = tic;
M0 = M(1:N,1:N);
toc(t0);
t0 = tic;
resp0 = M0 * exc;
toc(t0);

eps = log10(abs(resp./resp0-1));

figure;
plot_mesh(m, eps);
colorbar;
