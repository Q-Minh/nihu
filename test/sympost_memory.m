%% SYMPOST_MEMORY
% This script compares three algorithms for the solution of a symmetric
% postprocessing problems.
clear;
clc;

%% parameters
R = 1; % radius of sphere
Z = 1.2; % distance of sphere center from symmetry plane
L = 5; % receiver plane dimension
H = .2; % receiver plane height
D = 1; % receiver plane distance from sphere center
symm = 1; % rigid symmetry plane
ratio = 7; % bem element per wavelength ratio
shift = .3; % point source shifted from sphere origin
k = 20;

%% Create radiator and field point mesh
% radiator surface
Le = 2*pi/k/ratio; % element length
sphere = create_sphere_boundary(R, round(R/Le));
sphere = move_model(sphere, [0 0 Z]);
% field point mesh
field = create_slab(L*[1 1], round(L*[1 1]/Le));
field = rotate_model(field, [0 0 0], [1 0 0], pi/2);
field = move_model(field, [-mean(field.Nodes(:,2:3),1) H]);

%% Compute incident wave field on the bem mesh
r0 = [0 0 Z] + shift;
[source, trash, norm] = geo2gauss(sphere, [1 1]); % element centres
[pinc, qinc] = incident('point', r0, source, norm, k); % incident wave field

%% Compuatation with Different algorithms
% naive
[tree_naive, intdata_naive] = naive(sphere, field, k, [], [], symm);
% liu
[tree_liu, intdata_liu] = liu(sphere, field, k, [], [], symm);
% optimal
[tree_opti, intdata_opti] = opti(sphere, field, k, [], [], symm);

%% 
[B_naive, C_naive, M_naive] = fm_space(tree_naive, intdata_naive);
[B_liu, C_liu, M_liu] = fm_space(tree_liu, intdata_liu);
[B_opti, C_opti, M_opti] = fm_space(tree_opti, intdata_opti);

%%
f(1) = figure;
plot(0:length(tree_naive)-1, [B_naive, C_naive, M_naive]);
f(2) = figure;
plot(0:length(tree_liu)-1, [B_liu, C_liu, M_liu]);
f(3) = figure;
plot(0:length(tree_opti)-1, [B_opti, C_opti, M_opti]);
same_size(f);