function brick = create_brick_boundary(L, N)
%CREATE_BRICK_BOUNDARY  Create a brick surface mesh
%   BRICK  = CREATE_BRICK_BOUNDARY(L, N) creates a brick surface model with
%   side lengths given in L and division given in N. The brick is located
%   at the origin of the coordinate system, its faces are aligned along the
%   coordinate axes. The surface elements are aligned outward.
% Parameters:
%   L : 1x3 vector of side lengths of the brick mesh. If a scalar is given,
%       a cube model is created.
%   N : 1x3 vector of number of elements along the three directions. If a
%       scalar value is given, the number of elements is the same in the
%       three coordinate directions.
%
% See also: CREATE_BRICK, CREATE_LINE, CREATE_SLAB, CREATE_CIRCLE,
% CREATE_CIRCLE_QUADRANT, CREATE_SPHERE, CREATE_SPHERE_BOUNDARY,
% CREATE_CATSEYE

%   Copyright 2008-2010 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modifed: 02.12.2009

%% Parameter check
if nargin == 1
    N = L;
end
if isscalar(L)
    L = repmat(L,1,3);
end
if isscalar(N)
    N = repmat(N,1,3);
end

%% Create six faces of the brick model
% Create faces z = 0 and z = L(3)
s1 = create_slab(L(1:2), N(1:2));
s1 = join_meshes(flip_elements(s1), translate_mesh(s1, [0 0 L(3)]));
% Create faces x = 0 and x = L(1)
s2 = create_slab(L(3:-1:2), N(3:-1:2));
s2 = rotate_mesh(s2, [0 0 0], [0 1 0], -pi/2);
s2 = join_meshes(s2, translate_mesh(flip_elements(s2), [L(1) 0 0]));
% Create faces y = 0 and y = L(2)
s3 = create_slab(L([1 3]), N([1 3]));
s3 = rotate_mesh(s3, [0 0 0], [1 0 0], pi/2);
s3 = join_meshes(s3, translate_mesh(flip_elements(s3), [0 L(2) 0]));

%% Join faces
brick = join_meshes(s1, s2, s3);
brick = merge_coincident_nodes(brick);
brick = drop_unused_nodes(brick);
brick.Elements = drop_IDs(brick);
brick.Nodes(:,1) = 1:size(brick.Nodes,1);
end
