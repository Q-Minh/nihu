clear;

%%
profile on

%%
dim = 3;

if dim == 3
    msh = create_sphere_boundary(1, 100);
    depth = 7;
elseif dim == 2
    msh = create_circle(1, 30);
end
xc = centnorm(msh);
xc = xc(:,1:dim);
[T, fathersou, dsrc] = build_regular_clustertree(depth, xc);

%%
kernel = @laplace_kernel;
nExp = 4;
M2M = bb_M2M_regular(dim, nExp);
M2L = bb_M2L_regular(dim, kernel, nExp);
L2L = bb_L2L_regular(dim, nExp);
P2M = bb_P2M_regular(nExp, size(T(end).coord,1), fathersou, dsrc);
[i, j] = nfij(T(end).nearfield, T(end).nodsou, T(end).nodrec);
z = laplace_kernel_sp(xc(i,:), xc(j,:));
P2P = sparse(i, j, z);

%%
sigma = zeros(size(xc,1),1);
sigma(1) = 1;
sigma(round(end/2)) = -1;
sigma(round(end)) = 1;

tic;
fi_near = P2P * sigma;
[fi_far, times] = bb_far_transfer(T, sigma, P2M, M2M, M2L, L2L, P2M', -1);
fi = fi_near + fi_far;
toc;

%%
profile viewer

%%
figure;
plot_mesh(msh, 20*log10(fi));
h = findall(gca, 'type', 'patch');
set(h, 'linestyle', 'none');
set(gcf, 'renderer', 'painters');