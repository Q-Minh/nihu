clear;

%%
profile on

dim = 3;

if dim == 3
    msh = create_sphere_boundary(1, 50);
    xc = centnorm(msh);
    xc = xc(:,1:dim);
    [T, fathersou, dsrc] = build_regular_clustertree(6, xc);
elseif dim == 2
    msh = create_circle(1, 30);
end

%%
kernel = @laplace_kernel;
nExp = 4;
M2M = bb_M2M_regular(dim, nExp);
M2L = bb_M2L_regular(dim, kernel, nExp);
L2L = bb_L2L_regular(dim, nExp);
P2M = bb_P2M_regular(nExp, size(T(end).coord,1), fathersou, dsrc);

%%
sigma = zeros(size(xc,1),1);
sigma(1) = 1;
sigma(round(end/2)) = -1;
sigma(round(end)) = 1;

[fi, times] = bb_far_transfer(T, sigma, P2M, M2M, M2L, L2L, P2M', -1);

%%
figure;
subplot(1,2,1);
plot_mesh(msh, 20*log10(sigma));
h = findall(gca, 'type', 'patch');
set(h, 'linestyle', 'none');
subplot(1,2,2);
plot_mesh(msh, 20*log10(fi));
h = findall(gca, 'type', 'patch');
set(h, 'linestyle', 'none');
set(gcf, 'renderer', 'painters');