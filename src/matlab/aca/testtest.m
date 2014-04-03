clear;

m = create_sphere_boundary(1, 20);
k = min(mesh_kmax(m));
x = m.Nodes(:,2:4);

tic;
Ctree = build_cluster_tree(x);
tCtree = toc;
fprintf('%g s needed to build cluster tree\n', tCtree);

tic;
Btree = build_block_tree(Ctree);
tBtree = toc;
fprintf('%g s needed to build block tree\n', tBtree);

%%
N = size(x,1);

%%
q = 0;
qq = 0;
for b = 1 : length(Btree)
    bt = Btree(b,:);
    s = [length(Ctree(bt(1)).ind) length(Ctree(bt(2)).ind)];
    if any(s <= 2)
        q = q + 1;
        qq = qq + prod(s);
    end
end
fprintf(1, '%9d / %9d (%g) sparse entries\n', qq, N^2, qq/N^2);
fprintf(1, '%9d / %9d nontrivial blocks\n', length(Btree)-q, length(Btree));

%%
exc = ones(N,1);
resp = zeros(N,1);

M = @(row,col)helmholtz_matrix(row, col, x, k);

tic;
R = 3;
for b = 1 : size(Btree,1)
    ACA(b).I = Ctree(Btree(b,1)).ind;
    ACA(b).J = Ctree(Btree(b,2)).ind;
    [ACA(b).U, ACA(b).V] = lowrank_approx(M, ACA(b).I, ACA(b).J, 1e-3);
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
M0 = M(1:N,1:N);

eps = log10(abs(resp./resp0-1));

%%
figure;
plot_mesh(m, eps);
colorbar;
