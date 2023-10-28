function sphere = create_sphere(R, nR)
%CREATE_SPHERE Create a sphere volume mesh (NiHu / meshing)
%   SPHERE  = CREATE_SPHERE(R, N) creates a sphere volume model with radius
%   R and division N. The sphere is created as a "blown up cube" with
%   parameters L = [2 2 2] and Ncube = [2N 2N 2N], so N means the number of
%   elements along the radius
%
% See also: CREATE_LINE, CREATE_SLAB, CREATE_CIRCLE, CREATE_CIRCLE_QUADRANT,
% CREATE_BRICK, CREATE_BRICK_BOUNDARY, CREATE_SPHERE_BOUNDARY,
% CREATE_CATSEYE

%   Copyright 2008-2012 P. Fiala, P. Rucz
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modifed: 2012.12.14.

sphere = create_brick_base(nR*2*[1 1 1]);
nodes = sphere.Nodes(:,2:4);
nodes = bsxfun(@times, nodes, max(abs(nodes), [], 2)./sqrt(dot(nodes,nodes,2)));
nodes(isnan(nodes)) = 0;
sphere.Nodes(:,2:4) = nodes * R;

end
