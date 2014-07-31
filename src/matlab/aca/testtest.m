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
m = create_sphere_boundary(1, 24);
x = centnorm(m);

%%
t0 = tic;
Ctree = build_cluster_tree(x, 50);
tCtree = toc(t0);
fprintf('%.3g s needed to build cluster tree\n', tCtree);

%% plot cluster tree
v = ones(size(x,1),1);
b = 1;
for c = 1 : length(Ctree)
    if Ctree(c).level == 9
        v(Ctree(c).ind) = mod(b-1,20)+2;
        b = b+1;
    end
end
r = randperm(21)';
plot_mesh(m, r(v));

%% build block tree
t0 = tic;
[B_near, B_far] = build_block_tree(Ctree);
tBtree = toc(t0);
fprintf('%.3g s needed to build block tree\n', tBtree);

%%
figure;
display_block_structure(Ctree, B_near);

figure;
display_block_structure(Ctree, B_far);

% tic;
% M_near = sparse(S(:,1), S(:,2), Msp(S(:,1), S(:,2)));
% tNear = toc;
% fprintf('%.3g s needed to compute near field matrix\n', tNear);

%%
k = min(mesh_kmax(m));
M = @(i,j)helmholtz_matrix(i,j,x,k);

N = size(x,1);
exc = ones(N,1);

resp = zeros(N,1);

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
M0 = M(1:N,1:N);
toc(t0);
t0 = tic;
resp0 = M0 * exc;
toc(t0);

eps = log10(abs(resp./resp0-1));

figure;
plot_mesh(m, eps);
colorbar;
