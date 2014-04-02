clear;

m = create_sphere_boundary(1, 2);
x = m.Nodes(:,2:4);
tic;
Ctree = build_cluster_tree(x);
tCtree = toc;
fprintf('%g s needed to build cluster tree\n', tCtree);
tic;
Btree = build_block_tree(Ctree);
tBtree = toc;
fprintf('%g s needed to build block tree\n', tBtree);

N = size(x,1);
exc = ones(N,1);
resp = zeros(N,1);

M = @(row,col)laplace_matrix(row, col, x);

tic;
R = 3;
for b = 1 : length(Btree)
    ACA(b).I = Ctree(Btree(b,1)).ind;
    ACA(b).J = Ctree(Btree(b,2)).ind;
    [ACA(b).U, ACA(b).V] = lowrank_approx(M, ACA(b).I, ACA(b).J, R);
%     if any(any(isnan(ACA(b).U))) || any(any(isnan(ACA(b).V)))
%         aa = 2;
%     end
end
tACA = toc;
fprintf('%g s needed to generate ACA\n', tACA);

tic;
resp(:) = 0;
for b = 1 : length(Btree)
    I = ACA(b).I;
    J = ACA(b).J;
    resp(I) = resp(I) + ACA(b).U * (ACA(b).V * exc(J));
end
tmat = toc;
fprintf('%g s needed to compute matrix-vector product\n', tmat);

resp0 = M(1:N, 1:N) * exc;

eps = log10(abs(resp./resp0-1));

plot_mesh(m, eps);
colorbar;
% C = caxis(gca);
% C(1) = max(C(1), -6);
% caxis(gca, C);
% set(gca, 'clim', [-5 -1]);
