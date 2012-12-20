%% Exterior Neumann radiation problem
% This tutorial shows how to use the toolbox |NIHU| to compute a simple
% acoustic radiation problem with Neumann boundary condition.

%% Mesh generation
% Use a toolbox function |create_sphere_boundary| to create a sphere mesh
% with given radius and division parameter. Determine maximal wavelength
% applicable to the mesh with a given element-per-wavelength ratio, using
% the toolbox function |bemkmax|
R = 1;      % radius
nR = 20;    % division parameter
mesh = create_sphere_boundary(R, nR);
ratio = 7;  % element-per-wavelength ratio
kmax = mesh_kmax(mesh, ratio);
k = min(kmax);

figure;
plot_mesh(mesh, kmax); view(3);

%% Excitation
% The excitation will be the pressure derivative (velocity) field of a
% point source located inside the sphere, slightly shifed from the center.
% Constant elements are used with one DOF located at the element centers,
% so the velocity is computed in the element center locations. For error
% analysis purposes, the analytical pressure values are also
% computed.
% The element centers and normals are determined by
% a toolbox function |centnorm|, and the incident wave field is computed by
% the toolbox function |incident|.
r0 = [.2 0 0];                      % location of the source monopole
[cent, normal] = centnorm(mesh);   % element centers and normals
[ps0, qs] = incident('point', r0, cent, normal, k);

figure;
plot_mesh(mesh, qs); view(3); colorbar;

%% Build cluster tree
% The cluster tree is built with the toolbox function |clustertree|. The
% depth of the tree is determined by the mesh's dimensions and
% the minimal dimensionless cluster size at the leaf level as
%
% $$ depth = \log_2 \frac{2kR}{kd_{\mathrm{min}}}$$
%
kdmin = 3;  % minimal dimensionless cluster size at leaf level
kD = 2*k*R; % mesh dimension
depth = round(log2(kD/kdmin));
[tree, fs, fr] = clustertree(depth, cent);

%% Sparse matrices
% The near field sparse bem matrices are determined by the toolbox function
% |bemHG|. The sparsity structure (location of nonzeros) of these matrices
% is determined using the cluster tree by the toolbox function |nfij|
[i, j] = nfij(tree(end).nearfield, tree(end).nodsou);
[Hnf, Gnf] = bemHG(mesh, k, 'const', [], [i j]);

figure;
spy(Hnf);

%% FMM integration parameters
% The expansion lengths, the quadrature points over the unit sphere and the
% precomputed translation operators and interpolation matrices are
% generated by the toolbox function |integpar|.
C = 3;  % Accuracy parameter related to the expansion length
I = integpar(tree, k, C);

%% BEM spatial quadrature and relative tree
% The Gaussian integration points and weights on the surface are computed
% by the toolbox function |geo2gauss|. The parameters of the function are
% the mesh and the quadrature points for triangles and quadrangles,
% respectively.
% In order to optimally use the C (MEX) code, the tree is converted into a
% relative tree that contains distances between entries rather than
% absolute locations.
[gs, gn, w, gind] = geo2gauss(mesh, 3);
[tr, rr, rs, ns] = reltree(tree, cent, fr, gs, fs(gind), gn);

%% Solution
% The right hand side of the linear system of equations is computed now.
% The near field part is computed by multiplication by the near field
% matrix, the far field contribution is computed by means of the MLFMA,
% using the toolbox function |mpcont_Gq|.
Gq = Gnf * qs + ...
    mpcont_Gq(rr, rs, qs(gind).*w, tr, I, k, fs(gind), fr);

%% 
% The iterative solution of the system of equations is computed by the
% toolbox function |my_gmres|, where the matrix-vector products are
% evaluated by the toolbox function |gmres_iter_neumann|.
m = 50;     % max. number of iterations
tol = 1e-3; % prescribed backward relative tolerance
[ps, eps] = my_gmres(...
    @(p) gmres_iter_neumann(p, Hnf, rr, rs, ns, w, gind, tr, I, k, fs, fr),...
    Gq, m, tol, speye(size(Hnf)));

%%
% The solution is plotted with the toolbox function |plot_mesh|, and is
% compared with the analytical solution
e = norm(ps-ps0)/norm(ps0); % relative error

figure;
plot_mesh(mesh, real(ps)); view(3);
title(sprintf('Error: %.2e', e));
