function mesh = translate_mesh(mesh, dir)
%TRANSLATE_MESH Translate mesh along a given vector (NiHu / meshing)
%   MESH = TRANSLATE_MESH(MESH, DIR) translates the mesh MESH along the
%   direction vector DIR.
%
% See also: SCALE_MESH, ROTATE_MESH, REVOLVE_MESH, EXTRUDE_MESH,
% REPEAT_MESH, REFLECT_MESH

%   Copyright 2008-2010 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modifed: 2012.12.14.

ind = 1+(1:length(dir));
mesh.Nodes(:,ind) = bsxfun(@plus, mesh.Nodes(:,ind), dir(:).');

end
