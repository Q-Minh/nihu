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
D = 1; % receiver plane distance from sphere center
symm = 1; % rigid symmetry plane
ratio = 7; % bem element per wavelength ratio
shift = .3; % point source shifted from sphere origin
k = 40;

%% Create radiator and field point mesh
% radiator surface
Le = 2*pi/k/ratio; % element length
sphere = create_sphere_boundary(R, round(R/Le));
sphere = move_model(sphere, [0 0 Z]);
% field point mesh
field = create_slab(L*[1 1], round(L*[1 1]/Le));
field = rotate_model(field, [0 0 0], [1 0 0], pi/2);
field = move_model(field, [-mean(field.Nodes(:,2:3),1) H]);
[nodind, elind] = fe_select(field, sprintf('x.^2+(z-%g).^2 > %g', Z, (1.1*D)^2), 'ind');
field.Elements = field.Elements(elind,:);
drop_unused_nodes(field);

%% Compute incident wave field on the bem mesh
r0 = [0 0 Z] + shift;
[source, trash, norm] = geo2gauss(sphere, [1 1]); % element centres
[pinc, qinc] = incident('point', r0, source, norm, k); % incident wave field

%% Compuatation with Different algorithms
% naive
[tree_naive, intpar_naive, pf_naive, t_naive] = naive(sphere, field, k, pinc, qinc, symm);
% liu
[tree_liu, intpar_liu, pf_liu, t_liu] = liu(sphere, field, k, pinc, qinc, symm);
% optimal
[tree_opti, intpar_opti, pf_opti, t_opti] = opti(sphere, field, k, pinc, qinc, symm);
% analytical
pf_anal = incident('point', r0, field.Nodes(:,2:4), [], k) +...
    symm * incident('point', r0*diag([1 1 -1]), field.Nodes(:,2:4), [], k);

%% Compare different solutions
figure;
subplot(2,2,1);
plot_fem_model(field, field.Nodes(:,1), (abs(pf_naive)));
view(0,0);
subplot(2,2,2);
plot_fem_model(field, field.Nodes(:,1), (abs(pf_liu)));
view(0,0);
subplot(2,2,3);
plot_fem_model(field, field.Nodes(:,1), (abs(pf_opti)));
view(0,0);
subplot(2,2,4);
plot_fem_model(field, field.Nodes(:,1), (abs(pf_anal)));
view(0,0);

%% Compare different solutions
pinc = incident('point', r0, source, norm, k) +...
    symm * incident('point', r0*diag([1 1 -1]), source, norm, k);

figure;
p = patch('Vertices', sphere.Nodes(:,2:4), 'Faces', sphere.Elements(:,5:8),...
    'FaceVertexCData', abs(pinc));
shading faceted;
plot_fem_model(field, field.Nodes(:,1), (abs(pf_opti)));
shading faceted;
hold on;
axis equal;
axis tight;
view(-150, 20);
print -depsc figura.eps
