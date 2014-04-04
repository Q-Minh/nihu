clear;

m = create_sphere_boundary(1, 5);
k = min(mesh_kmax(m));
x = m.Nodes(:,2:4);
M = @(row,col)helmholtz_matrix(row, col, x, k);
Msp = @(row,col)helmholtz_matrix_sp(row, col, x, k);

tic;
Ctree = build_cluster_tree(x);
tCtree = toc;
fprintf('%.3g s needed to build cluster tree\n', tCtree);

tic;
[B_near, B_far] = build_block_tree(Ctree);
tBtree = toc;
fprintf('%.3g s needed to build block tree\n', tBtree);

% tic;
% M_near = sparse(S(:,1), S(:,2), Msp(S(:,1), S(:,2)));
% tNear = toc;
% fprintf('%.3g s needed to compute near field matrix\n', tNear);

tic;
for b = 1 : size(B_far,1)
    ACA(b).I = Ctree(B_far(b,1)).ind;
    ACA(b).J = Ctree(B_far(b,2)).ind;
    [ACA(b).U, ACA(b).V] = lowrank_approx(M, ACA(b).I, ACA(b).J, 1e-3);
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
