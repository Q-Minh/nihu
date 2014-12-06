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
Ns = 5000;
Nr = 7000;
source = create_line(1, Ns);
receiver = translate_mesh(create_line(1, Nr), [.5 0 0]);
xs = centnorm(source);
xr = centnorm(receiver);

%% build cluster trees
t0 = tic;
nLeaf = 25;
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
for b = 1 : size(B_far,1)
    i = ReceiverClusters(B_far(b,1),:);
    j = SourceClusters(B_far(b,2),:);
    patch(j(1) + j(2)*[0 0 1 1]', i(1) + i(2)*[0 1 1 0]', b*ones(4,1));
end

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
k = min(mesh_kmax(source));
M = @(i,j)helmholtz_matrix(i,j,xr,xs,k);

exc = ones(Ns,1);
resp = zeros(Nr, 1);

t0 = tic;
for i = 1 : size(B_near,1)
    I = ReceiverTree(B_near(i,1)).ind;
    J = SourceTree(B_near(i,2)).ind;
    resp(I) = resp(I) + M(I,J) * exc(J);
end
tNear = toc(t0);
fprintf('%.3g s needed to compute near field\n', tNear);

eps = 1e-2;
t0 = tic;
for tau = 1 : size(B_far,1)
    progbar(1, size(B_far,1), tau);
    I = ReceiverTree(B_far(tau,1)).ind;
    J = SourceTree(B_far(tau,2)).ind;
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
plot(xr(:,1), logeps);
xlabel('x [m]');
ylabel('log10 error [-]');
ylim([-10 0]);
