function model = create_circle_quadrant(R, nR)
%CREATE_CIRCLE_QUADRANT Create quadrant of a circle surface mesh
%   CIRCLE = CREATE_CIRCLE_QUADRANT(R, N) creates the quadrant of a
%   circular mesh located at the origin. The radius of the circle is R, its
%   division parameter is N. N is equal to the number of elements along the
%   radius, this value needs to be even, otherwise it is modified to the
%   nearest larger even number.
%   The circle quadrant is located in the xy plane, its sides are parallel
%   with the x and y axes.
%
% See also: CREATE_LINE, CREATE_SLAB, CREATE_CIRCLE, CREATE_BRICK,
% CREATE_BRICK_BOUNDARY, CREATE_SPHERE, CREATE_SPHERE_BOUNDARY,
% CREATE_CATSEYE

%   Copyright 2008-2012 P. Fiala, P. Rucz
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modifed: 2012.12.14.


nR = ceil(nR/2)*2; % number of elements over the radius (have to be even)
center = .6/sqrt(2) * [1 1]; % 'center' point of the quadrant
% internal slice
cin = [
     0 0
    .5 0
    center
     0 .5
    ];
inner = create_slab(cin, nR/2*[1 1]);
% outer slice
outer = create_slab([1 1], nR/2*[1 1]);
level = outer.Nodes(:,2);
cou = [
    .5 0
    1 0
    1/sqrt(2) 1/sqrt(2)
    center
    ];
outer1 = create_slab(cou, nR/2*[1 1]);
[cou(:,2), cou(:,1)] = cart2pol(cou(:,1), cou(:,2));
outer2 = create_slab(cou, nR/2*[1 1]);
[outer2.Nodes(:,2), outer2.Nodes(:,3)] = pol2cart(outer2.Nodes(:,3), outer2.Nodes(:,2));
level = repmat(level,1,3);
outer.Nodes(:,2:4) = (1-level) .* outer1.Nodes(:,2:4) + level .* outer2.Nodes(:,2:4);
% combine innner, outer and outer's image regions
model = join_meshes(inner, outer, reflect_mesh(outer, [-1 +1 0]));
% scale to radius
model.Nodes(:,2:4) = R * model.Nodes(:,2:4);
% merge coincident nodes
tol = R / (10*nR);
model = merge_coincident_nodes(model, tol);
model = drop_mesh_IDs(model);

end
