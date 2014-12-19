clear;

%%
profile on

%%
dim = 2;

if dim == 3
    msh = create_sphere_boundary(1, 40);
    depth = 6;
elseif dim == 2
    msh = create_circle(1, 70);
    depth = 6;
end
xc = centnorm(msh);
xc = xc(:,1:dim);
[T, fathersou, dsrc] = build_regular_clustertree(depth, xc);

%%
kernel = @laplace_kernel;
kernel_sp = @laplace_kernel_sp;
nExp = 5;
P2P = bb_P2P_regular(T(end), xc, xc, kernel_sp);
P2M = bb_P2M_regular(nExp, size(T(end).coord,1), fathersou, dsrc);
M2M = bb_M2M_regular(dim, nExp);
M2L = bb_M2L_regular(dim, kernel, nExp);
L2L = bb_L2L_regular(dim, nExp);

% eps = 1e-4;
% [CM2L.U, CM2L.V, CM2L.r] = compress_M2L(M2L, eps);

%%
matfun = @(x)bb_matvec_regular(x, T, P2P, P2M, M2M, M2L, L2L, P2M.', -1);
L = .2;
sigma0 = xc(:,1);
fi0 = matfun(sigma0);
[sigma,flag,relres,iter,resvec] = gmres(matfun, fi0);

%%
profile viewer

%%
figure;
plot_mesh(msh, sigma);
h = findall(gca, 'type', 'patch');
set(h, 'linestyle', 'none');
set(gcf, 'renderer', 'painters');
