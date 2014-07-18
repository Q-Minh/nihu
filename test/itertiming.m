%% SYMPOST
% This script compares three algorithms for the solution of a symmetric
% postprocessing problems.
clear;
clc;

%% parameters
R = 1; % radius of sphere
Z = 1.2; % distance of sphere center from symmetry plane
L = 5; % receiver plane dimension
H = .2; % receiver plane height
symm = 1; % rigid symmetry plane
ratio = 7; % bem element per wavelength ratio
shift = .3; % point source shifted from sphere origin
k = 60;

%% Create radiator and field point mesh
% radiator surface
Le = 2*pi/k/ratio; % element length
model = create_sphere_boundary(R, round(R/Le));
model = move_model(model, [0 0 Z]);

%% Compute incident wave field on the bem mesh
r0 = [0 0 Z] + shift;
[ecoord, trash, normal] = geo2gauss(model, [1 1]); % element centres
[pinc, qinc] = incident('point', r0, ecoord, normal, k); % incident wave field

%% Gaussian integration points normals and weights
ng = [0 16];
[ecoord, gcoord, gnorm, w, gind] = geo2gauss(model, ng); % Quadrature config

%% cluster tree
% determine tree depth
D = norm(max(model.Nodes(:,2:4),[],1) - min(model.Nodes(:,2:4),[],1));
depth = round(log2(k*D/pi));
% build tree
ttstart = tic;
[tree fathersou fatherrec] = clustertree2(depth, ecoord, [], symm);  % build cluster tree using elem centers
print_tree_info(tree);
% Gaussian point fathers
fathergau = zeros(size(gcoord,1),1);
fathergau(gind) = repmat(fathersou,1,size(gind,2));
tt = toc(ttstart);

%% FMBEM integration parameters
C = 4;
intdata = integpar(tree, k, C, symm);     % determine integration parameters
init_translation(tree, intdata, k, 0); % store translation operators

%% Right hand side
% Excitation in Gaussian points
pw = repmat(pinc.',length(w)/length(pinc),1);
pw = pw(:) .* w;
% multipole contribution
[Hpm tleaf ttrans tsampup tshiftup tsampdn tshiftdn trec] = mpcont_Hp2(ecoord, gcoord, gnorm, pw, tree, intdata, k, fathergau, fatherrec, symm);
%%
t = [tshiftup tsampup ttrans tshiftdn tsampdn];
bar(t, 'Stack');
legend({'Shift up', 'Interpolate', 'Translation', 'Shift down', 'Filter'}, 'Location', 'NorthWest');
