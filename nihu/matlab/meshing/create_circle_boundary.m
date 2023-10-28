function model = create_circle_boundary(R, n)
%CREATE_CIRCLE_BOUNDARY Create a circle boundary mesh (NiHu / meshing)
%   CIRCLE = CREATE_CIRCLE_BOUNDARY(R, N) creates a circle mesh located
%   at the origin. The radius of the circle is R, and the circle is
%   divided into N segments. The circle is located in the xy plane.
%
% See also: CREATE_LINE, CREATE_SLAB, CREATE_CIRCLE_QUADRANT, CREATE_BRICK,
% CREATE_BRICK_BOUNDARY, CREATE_SPHERE, CREATE_SPHERE_BOUNDARY,
% CREATE_CATSEYE

%   Copyright 2008-2013 P. Fiala, P. Rucz
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modified 2013.08.14.

% Create and reflect a circle quadrant
model = create_line(2*pi, n);
phi = model.Nodes(:,2);
model.Nodes(:,2:3) = R * [cos(phi) sin(phi)];

% Merge model and get rid of IDs
model = merge_coincident_nodes(model, R/n/10);
model = drop_mesh_IDs(model);

end
