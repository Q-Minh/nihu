clear;

%%
dim = 2;

if dim == 3
    msh = create_sphere_boundary(1, 30);
elseif dim == 2
    msh = create_circle(1, 30);
end
xc = centnorm(msh);
xc = xc(:,1:dim);
T = build_regular_clustertree(4, xc);

%%
kernel = @laplace_kernel;
nExp = 4;
M2M = bb_M2M_regular(dim, nExp);
M2L = bb_M2L_regular(dim, kernel, nExp);
L2L = bb_L2L_regular(dim, nExp);

%%
sigma = ones(nExp^dim, size(T(end).coord,1));
bb_far_transfer(T, sigma, M2M, M2L, L2L);
