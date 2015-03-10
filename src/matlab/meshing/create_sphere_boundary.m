function mesh = create_sphere_boundary(R, nR)
%CREATE_SPHERE_BOUNDARY Create a sphere surface mesh (NiHu / meshing)
%   SPHERE  = CREATE_SPHERE_BOUNDARY(R, N) creates a sphere surface mesh
%   with radius R and division N. The sphere is created as a "blown up
%   cube surface" with parameters L = [2 2 2] and Ncube = [2N 2N 2N], so N
%   means the number of elements along the radius.
%
% See also: CREATE_LINE, CREATE_SLAB, CREATE_CIRCLE, CREATE_CIRCLE_QUADRANT,
% CREATE_BRICK, CREATE_BRICK_BOUNDARY, CREATE_SPHERE, CREATE_CATSEYE

%   Copyright 2008-2012 P. Fiala, P. Rucz
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modifed: 2012.12.14.


mesh = create_brick_boundary_base(nR*2*[1 1 1]);
nodes = mesh.Nodes(:,2:4);
nodes = bsxfun(@times, nodes, max(abs(nodes), [], 2)./sqrt(dot(nodes,nodes,2)));
mesh.Nodes(:,2:4) = nodes * R;

end
