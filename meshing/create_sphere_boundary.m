function mesh = create_sphere_boundary(R, nR)
%CREATE_SPHERE_BOUNDARY  Create a sphere surface mesh
%   SPHERE  = CREATE_SPHERE_BOUNDARY(R, N) creates a sphere surface mesh
%   with radius R and division N. The sphere is created as a "blown up
%   cube surface" with parameters L = [2 2 2] and Ncube = [2N 2N 2N], so N
%   means the number of elements along the radius.
%
% See also: CREATE_LINE, CREATE_SLAB, CREATE_CIRCLE, CREATE_CIRCLE_QUADRANT,
% CREATE_BRICK, CREATE_BRICK_BOUNDARY, CREATE_SPHERE, CREATE_CATSEYE

%   Copyright 2008-2010 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modifed: 02.12.2009

%% Parameter check
error(nargchk(2, 2, nargin, 'struct'));

%%
mesh = translate_mesh(create_brick_boundary(2, nR*2), -[1 1 1]);
nodes = mesh.Nodes(:,2:4);
for k = 1 : 3
    nodes(:,k) = tan(pi/4*nodes(:,k));
end
nodes = nodes ./ repmat(sqrt(dot(nodes,nodes,2)),1,3);
mesh.Nodes(:,2:4) = R * nodes;
end
