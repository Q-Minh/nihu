function mesh = reflect_mesh(mesh, normal, base)
%REFLECT_MESH Reflect mesh to a symmetry plane (NiHu / meshing)
%   MESH = REFLECT_MESH(MESH, NORMAL) or 
%   MESH = REFLECT_MESH(MESH, NORMAL, BASE) reflects the NiHu mesh MESH to
%   the symmetry plane defined by BASE and NORMAL.
%   NORMAL : normal of the symmetry plane
%   BASE   : a point of the symmetry plane, default is [0 0 0]
%
% See also: TRANSLATE_MESH, SCALE_MESH, ROTATE_MESH, EXTRUDE_MESH,
% REVOLVE_MESH, REPEAT_MESH

%   Copyright 2008-2012 P. Fiala, P. Rucz
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last updated: 2012.12.14.

% Argument check
if nargin < 3
    base = [0 0 0];
else
	base = base(:).';
end

normal = normal(:).' / norm(normal);

% reflect nodes
nodes = mesh.Nodes(:,2:4);
mesh.Nodes(:,2:4) = nodes - 2*(bsxfun(@minus, nodes, base) * normal.') * normal;

% flip elements
mesh.Elements = flip_elements(mesh.Elements);
end
