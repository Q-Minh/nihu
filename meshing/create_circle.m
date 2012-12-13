function model = create_circle(R, nR, varargin)
%CREATE_CIRCLE Create a circle surface mesh
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

% Last modified 2012.12.13.

% Create and reflect a circle quadrant
model = create_circle_quadrant(R, nR, varargin{:});
model = join_meshes(model, reflect_mesh(model, [0 0 0], [1 0 0]));
model = join_meshes(model, reflect_mesh(model, [0 0 0], [0 1 0]));

% Merge model and get rid of IDs
model = merge_coincident_nodes(model,varargin{:});
model.Elements = drop_IDs(model);
model.Elements(:,[3 4]) = repmat([1,1],size(model.Elements,1),1);
model.Properties = [1 1];
model.Materials = [1 1];
model.Nodes(:,1) = 1 : size(model.Nodes,1);
end
