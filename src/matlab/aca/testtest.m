clear;

%%
if isunix
    root = '/D';
else
    root = 'D:';
end
m = import_mesh(fullfile(root, 'research', 'pub',...
    '2013', 'Boonen2013', 'work', 'industrial', 'data',...
    'horse.off'));

%%
% m = create_sphere_boundary(1, 15);
x = centnorm(m);

%%
tic;
Ctree = build_cluster_tree(x);
tCtree = toc;
fprintf('%.3g s needed to build cluster tree\n', tCtree);

%% plot cluster tree
v = ones(size(x,1),1);
b = 1;
for c = 1 : length(Ctree)
    if Ctree(c).level == 13
        v(Ctree(c).ind) = mod(b-1,20)+2;
        b = b+1;
    end
end
r = randperm(21)';
plot_mesh(m, r(v));

%%
tic;
[B_near, B_far] = build_block_tree(Ctree);
tBtree = toc;
fprintf('%.3g s needed to build block tree\n', tBtree);

figure;
display_block_structure(Ctree, B_near);
figure;
display_block_structure(Ctree, B_far);

% tic;
% M_near = sparse(S(:,1), S(:,2), Msp(S(:,1), S(:,2)));
% tNear = toc;
% fprintf('%.3g s needed to compute near field matrix\n', tNear);

tic;
for b = 1 : size(B_far,1)
    ACA(b).I = Ctree(B_far(b,1)).ind;
    ACA(b).J = Ctree(B_far(b,2)).ind;
    [ACA(b).U, ACA(b).V] = lowrank_approx_block(M, ACA(b).I, ACA(b).J, 1e-3);
end
tACA = toc;
fprintf('%.3g s needed to generate ACA\n', tACA);

N = size(x,1);
exc = ones(N,1);

tic;
resp = M_near * exc;
for b = 1 : length(B_far)
    I = ACA(b).I;
    J = ACA(b).J;
    resp(I) = resp(I) + ACA(b).U * (ACA(b).V * exc(J));
end
tmat = toc;
fprintf('%.3g s needed to compute matrix-vector product\n', tmat);

tic;
M0 = M(1:N,1:N);
toc;
tic;
resp0 = M0 * exc;
toc;

eps = log10(abs(resp./resp0-1));

figure;
plot_mesh(m, eps);
colorbar;
