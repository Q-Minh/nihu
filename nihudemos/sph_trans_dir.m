%% Exterior Dirichlet radiation problem
% This tutorial shows how to use the toolbox |NIHU| to compute a simple
% acoustic radiation problem with Dirichlet boundary condition.

%% Mesh generation
% Use a toolbox function |create_sphere_boundary| to create a sphere mesh
% with given radius and division parameter. Determine maximal wavelength
% applicable to the mesh with a given element-per-wavelength ratio, using
% the toolbox function |bemkmax|
R = 1;      % radius
nR = 10;    % division parameter
mesh = create_sphere_boundary(R, nR);
ratio = 7;  % element-per-wavelength ratio
kmax = bemkmax(mesh, ratio);
k = min(kmax);

figure;
plot_mesh(mesh, kmax); view(3);

%% Excitation
% The excitation will be the pressure field of a point source
% located inside the sphere, slightly shifed from the center.
% Constant elements are used with one DOF located at the element centers,
% so the pressure is computed in the element center locations. For error
% analysis purposes, the analytical pressure derivative values are also
% computed.
% The element centers and normals are determined by
% a toolbox function |centnorm|, and the incident wave field is computed by
% the toolbox function |incident|.
r0 = [-.2 0 0];                    % location of the source monopole
[cent, normal] = centnorm(mesh);   % element centers and normals
[ps, qs0] = incident('point', r0, cent, normal, k);

figure;
plot_mesh(mesh, ps); view(3); colorbar;

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
%
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
Hp = Hnf*ps -ps/2 + ...
    mpcont_Hp(rr, rs, ns, ps(gind).*w, tr, I, k, fs(gind), fr);

%% 
% The iterative solution of the system of equations is computed by the
% toolbox function |my_gmres|, where the matrix-vector products are
% evaluated by the toolbox function |gmres_iter_dirichlet|.
m = 100;        % max. number of iterations
tol = 1e-3;     % prescribed backward relative tolerance
[qs, eps] = my_gmres(...
    @(q) gmres_iter_dirichlet(q, Gnf, rr, rs, w, gind, tr, I, k, fs, fr),...
    Hp, m, tol, speye(size(Gnf)));

%% 
% The solution is plotted with the toolbox function |plot_mesh|, and is
% compared with the analytical solution
e = norm(qs-qs0)/norm(qs0); % relative error

figure;
plot_mesh(mesh, real(qs)); view(3);
title(sprintf('Error: %.2e', e));

%% Field point mesh
xlimits = R * [1.1 3];
tlimits = [0 pi/2];
line = create_line([xlimits(1) 0 0; xlimits(2) 0 0], ceil(diff(xlimits)*k*ratio/2/pi));
ntheta = ceil(max(xlimits)*diff(tlimits)*k*ratio/2/pi);
dtheta = diff(tlimits)/ntheta;
field = revolve_mesh(line, [0 0 0], [0 -1 0], dtheta, ntheta);

%% Build cluster tree
points = field.Nodes(:,2:4);
[tree, fs, fr] = clustertree(depth, cent, points);

%% Near field
[i, j] = nfij(tree(end).nearfield, tree(end).nodsou, tree(end).nodrec);
[Hnf, Gnf] = bemHG(mesh, k, 'const', points, [i j]);

%%
C = 3;  % Accuracy parameter related to the expansion length
I = integpar(tree, k, C);

%%
[tr, rr, rs, ns] = reltree(tree, points, fr, gs, fs(gind), gn);

%%
pf = Hnf * ps + mpcont_Hp(rr, rs, ns, ps(gind).*w, tr, I, k, fs(gind), fr) -...
    Gnf * qs - mpcont_Gq(rr, rs, qs(gind).*w, tr, I, k, fs(gind), fr);

%% plot results
figure;
pf0 = incident('point', r0, points, [], k);
plot_mesh(field, pf);
shading interp;
plot_mesh(mesh, ps);
title(sprintf('Error: %.3x', norm(pf-pf0)/norm(pf0)));