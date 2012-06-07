function N = symmetric_fuso2(k)

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

N = [size(sphere.Nodes,1), size(field.Nodes,1)];

% %% compute incident wave field
% T = diag([1 1 -1]);
% [source, trash, norm] = geo2gauss(sphere, [1 1]);
% [pinc, qinc] = incident('point', r0, source, norm, k);
% [pinc2, qinc2] = incident('point', r0*T, source, norm, k);
% pinc = pinc + symm * pinc2;
% qinc = qinc + symm * qinc2;
% 
% %% Compute transfer with fmbem
% [pf t] = fm_postproc(sphere, field, k, qinc, pinc, symm, [3 4]);
% 
% % %% Analytical
% % points = field.Nodes(:,2:4);
% % panal = incident('point', r0, points, [], k) + ...
% %     symm * incident('point', r0*T, points, [], k);
% % 
% % %% plot results
% % figure;
% % subplot(2,1,1);
% % patch('Vertices', sphere.Nodes(:,2:4), 'Faces', sphere.Elements(:,5:8), ...
% %     'FaceVertexCData', abs(pinc));
% % shading faceted;
% % plot_fem_model(field, field.Nodes(:,1), log(abs(p)));
% % view(15,5);
% % axis equal; axis tight;
% % 
% % subplot(2,1,2);
% % patch('Vertices', sphere.Nodes(:,2:4), 'Faces', sphere.Elements(:,5:8), ...
% %     'FaceVertexCData', abs(pinc));
% % shading faceted;
% % plot_fem_model(field, field.Nodes(:,1), log(abs(panal)));
% % view(15,5);
% % axis equal; axis tight;