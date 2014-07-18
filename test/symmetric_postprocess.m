%% SYMMETRIC POSTPROCESS
% Short description
clear;
close;
clc;

%% Parameters
R = 1; % radius of source sphere
Z = 1.2; % distance from symmetry plane
D = 2; % distance of receiver plane from source
L = 5; % receiver plane side length
k = 40; % wave number
symm = 1; % positive symmetry
ratio = 7; % meshing ratio (d_elem / lambda)
shift = .3;
r0 = shift + [0 0 Z]; % source location

%% mesh
% radiator surface
Le = 2*pi/k/ratio;
singlesphere = translate_mesh(create_sphere_boundary(R, round(R/Le)), [0 0 Z]);
doublesphere = join_meshes(singlesphere, reflect_mesh(singlesphere, [0 0 0], [0 0 1]));

%% field point mesh
field = create_slab(L, round(L/Le));
field = translate_mesh(rotate_mesh(translate_mesh(field, [D 0 0]), [0 0 0], [1 0 0], pi/2), [0 0 Z]);
% reflect sphere to z = 0 plane

%% compute incident wave field
T = diag([1 1 -1]);
order = 3;
[cent, normal] = centnorm(singlesphere);
[pinc, qinc] = incident('point', r0, cent, normal, k);
[pinc2, qinc2] = incident('point', r0*T, cent, normal, k);
pinc = pinc + symm * pinc2;
qinc = qinc + symm * qinc2;
pinc2 = [pinc; symm*pinc];
qinc2 = [qinc; symm*qinc];

%% Compute transfer with symmetric and normal fmbem
[pf t] = fm_postproc(doublesphere, field, k, qinc2, pinc2, 0, order);
[pfs ts] = fm_postproc(singlesphere, field, k, qinc, pinc, symm, order);

%% Analytical
points = field.Nodes(:,2:4);
panal = incident('point', r0, points, [], k) + ...
    symm * incident('point', r0*T, points, [], k);

%% plot results
figure;
subplot(1,3,1);
plot_mesh(singlesphere, log(abs(pinc)));
plot_mesh(field, log(abs(pf)));

subplot(1,3,2);
plot_mesh(singlesphere, log(abs(pinc)));
plot_mesh(field, log(abs(pfs)));

subplot(1,3,3);
plot_mesh(singlesphere, log(abs(pinc)));
plot_mesh(field, log(abs(panal)));

setfig('linestyle', 'none');
