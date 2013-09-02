function mesh = create_catseye(R, nR)
%CREATE_CATSEYE Create a Cat's eye surface mesh (NiHu / meshing)
%   MESH = CREATE_CATSEYE(R, N) Creates a Cat's eye surface mesh with
%   radius R and division parameter N.
%
% See also: CREATE_BRICK_BOUNDARY, CREATE_LINE, CREATE_SLAB, CREATE_CIRCLE,
% CREATE_CIRCLE_QUADRANT, CREATE_SPHERE, CREATE_SPHERE_BOUNDARY,
% CREATE_CATSEYE

%   Copyright 2008-2010 P. Fiala, P. Rucz
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modifed: 2012.12.14.


% sphere boundary
mesh = create_sphere_boundary(R, nR);
% cut one octant
eps = 1e-2 * R/nR;
expression = sprintf('z < 0+%g | phi > pi/2-%g | phi < 0+%g', eps, eps, eps);
[nodind, elind] = mesh_select(mesh, expression, 'ind');
mesh.Elements = mesh.Elements(elind,:);
mesh.Nodes = mesh.Nodes(nodind,:);
% quadrant circle (xy)
quad = create_circle_quadrant(R, 2*nR);
% rotate twice to (xz)
quad2 = rotate_mesh(quad, [-pi/2; -pi/2]);
% rotate twice to (yz)
quad3 = rotate_mesh(quad, [pi/2; pi/2], [0 1 0; 1 0 0]);
% join meshs
mesh = join_meshes(mesh, quad, quad2, quad3);
% merge and drop IDs
mesh = merge_coincident_nodes(mesh, R/nR/10);
mesh = drop_mesh_IDs(mesh);

end
