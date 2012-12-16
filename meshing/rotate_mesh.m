function mesh = rotate_mesh(mesh, ang, dir, base)
%ROTATE_MESH Rotate mesh around a given vector (NiHu / meshing)
%   MESH = ROTATE_MESH(MESH, ANG, DIR, BASE) rotates the mesh MESH
%
% See also: TRANSLATE_MESH, SCALE_MESH, REVOLVE_MESH, EXTRUDE_MESH,
% REPEAT_MESH, REFLECT_MESH

%   Copyright 2008-2010 P. Fiala, P. Rucz
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modifed: 2012.12.14.

if nargin < 3
    dir = eye(numel(ang));
end
dir = [dir, zeros(size(dir,1), 3-size(dir,2))];

% construct rotation matrix
T = eye(3);
for i = 1 : numel(ang)
    T = T * rotation_matrix(ang(i), dir(i,:));
end

% perform rotation
if nargin > 3
    mesh = translate_mesh(mesh, -base);
end
mesh.Nodes(:,2:4) = mesh.Nodes(:,2:4) * T;
if nargin > 3
    mesh = translate_mesh(mesh, base);
end

end
