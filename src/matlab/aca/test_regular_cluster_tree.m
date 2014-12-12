clear;

%%
msh3 = create_sphere_boundary(1, 50);
xc3 = centnorm(msh3);
T3 = build_regular_clustertree(6, xc3);

dim = 3;
kernel = @laplace_kernel;
nExp = 4;
M2L = bb_M2L_regular(dim, kernel, nExp);

%%
msh2 = create_circle(1, 50);
xc2 = centnorm(msh2);
xc2 = xc2(:,1:2);
T2 = build_regular_clustertree(8, xc2);
