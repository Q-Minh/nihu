function t = symmetric_fuso(k)

%% Parameters
R = 1; % radius of source sphere
Z = 1.2; % distance from symmetry plane
D = 5; % distance of receiver plane from source
L = 5; % receiver plane side length
symm = 1; % positive symmetry
ratio = 7; % meshing ratio (d_elem / lambda)
shift = .3;
r0 = shift + [0 0 Z]; % source location

%% mesh
% radiator surface
Le = 2*pi/k/ratio;
sphere = create_sphere_boundary(R, round(R/Le));
sphere = move_model(sphere, [0 0 Z]);
% field point mesh
field = create_slab(L*[1 1], round(L*[1 1]/Le));
field = rotate_model(move_model(field, [D 0 0]), [0 0 0], [1 0 0], pi/2);

%% compute incident wave field
T = diag([1 1 -1]);
[source, trash, norm] = geo2gauss(sphere, [1 1]);
[pinc, qinc] = incident('point', r0, source, norm, k);
[pinc2, qinc2] = incident('point', r0*T, source, norm, k);
pinc = pinc + symm * pinc2;
qinc = qinc + symm * qinc2;

%% Compute transfer with fmbem
[pf t] = fm_postproc(sphere, field, k, qinc, pinc, symm, [3 4]);