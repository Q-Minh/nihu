function mesh = brick2cylinder(brick_mesh)
% BRICK2CYLINDER Transforms a brick mesh into a cylinder.
%   MESH = BRICK2CYLINDER(BRICK_MESH) returns a cylinder shaped mesh, which
%       is a transformation of the original mesh BRICK_MESH. The diameter
%       of the cylinder will be equal to the depth (or width) of the brick.

%   Copyright 2010 P. Rucz
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last updated: 06.09.2010.

mesh = brick_mesh;
xy = abs(mesh.Nodes(:,2:3));
phi = atan(min(xy,[],2)./max(xy,[],2));
scale = cos(phi);
scale(any(xy == 0, 2)) = 1;
mesh.Nodes(:,2:3) = mesh.Nodes(:,2:3).*([scale scale]);

