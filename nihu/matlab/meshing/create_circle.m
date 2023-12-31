function model = create_circle(R, nR, varargin)
%CREATE_CIRCLE Create a circle surface mesh (NiHu / meshing)
%   CIRCLE = CREATE_CIRCLE(R, N) creates a circular mesh located at the
%   origin. The radius of the circle is R, its division parameter is N. N
%   is equal to the number of elements along the radius, this value needs
%   to be even, otherwise it is modified to the nearest larger even number.
%   The circle is located in the xy plane.
%
% See also: CREATE_LINE, CREATE_SLAB, CREATE_CIRCLE_QUADRANT, CREATE_BRICK,
% CREATE_BRICK_BOUNDARY, CREATE_SPHERE, CREATE_SPHERE_BOUNDARY,
% CREATE_CATSEYE

%   Copyright 2008-2010 P. Fiala, P. Rucz
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modified 2012.12.14.

% Create and reflect a circle quadrant
model = create_circle_quadrant(R, nR);
model = join_meshes(model, reflect_mesh(model, [1 0 0]));
model = join_meshes(model, reflect_mesh(model, [0 1 0]));

% Merge model and get rid of IDs
model = merge_coincident_nodes(model, R/nR/10);
model = drop_mesh_IDs(model);

end
