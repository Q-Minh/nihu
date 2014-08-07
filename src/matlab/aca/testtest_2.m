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
source = create_sphere_boundary(1, 30);
receiver = translate_mesh(create_sphere_boundary(1, 30), [.5 0 0]);
xs = centnorm(source);
xr = centnorm(receiver);

%%
t0 = tic;
SourceTree = build_cluster_tree(xs, 25);
ReceiverTree = build_cluster_tree(xr, 25);
tTree = toc(t0);
fprintf('%.3g s needed to build cluster tree\n', tTree);

%%
l = 10;
ns = 0;
for s = 1: length(SourceTree)
    if SourceTree(s).level == l
        ns = ns + length(SourceTree(s).ind);
    end
end

nr = 0;
for r = 1: length(ReceiverTree)
    if ReceiverTree(r).level == l
        nr = nr + length(ReceiverTree(r).ind);
    end
end

%% build block tree
t0 = tic;
[B_near, B_far] = build_block_tree_2(SourceTree, ReceiverTree);
tBtree = toc(t0);
fprintf('%.3g s needed to build block tree\n', tBtree);

[~, i] = sort(B_near(:,1));
B_near = B_near(i,:);

[~, i] = sort(B_far(:,1));
B_far = B_far(i,:);

%% visualise far field block tree

v = nan(size(xr,1),1);
w = nan(size(xs,1),1);
sigma = unique(B_far(:,1));
sigma = sigma(end);
v(ReceiverTree(sigma).ind) = 1;
tau = B_far(B_far(:,1) == sigma,2);
for t = 1: length(tau);
    w(SourceTree(tau(t)).ind) = t;
end

figure;
plot_mesh(source, w);
plot_mesh(receiver, v);

shading flat;

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
