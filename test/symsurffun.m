function t = symsurffun(k)

%% parameters
R = 1; % radius of sphere
Z = 1.2; % distance of sphere center from symmetry plane
symm = 1; % rigid symmetry plane
ratio = 7; % bem element per wavelength ratio
shift = .3; % point source shifted from sphere origin

%% Create radiator and field point mesh
% radiator surface
Le = 2*pi/k/ratio; % element length
sphere = create_sphere_boundary(R, round(R/Le));
sphere = move_model(sphere, [0 0 Z]);

%% Compute incident wave field on the bem mesh
r0 = [0 0 Z] + shift;
[source, trash, norm] = geo2gauss(sphere, [1 1]); % element centres
[pinc, qinc] = incident('point', r0, source, norm, k); % incident wave field

%% Compuatation with Different algorithms
% naive
[trash, trash, trash, t(1,:)] = naive_surf(sphere, k, pinc, qinc, symm);
% liu
[trash, trash, trash, t(2,:)] = liu_surf(sphere, k, pinc, qinc, symm);
% optimal
[trash, trash, trash, t(3,:)] = opti_surf(sphere, k, pinc, qinc, symm);
