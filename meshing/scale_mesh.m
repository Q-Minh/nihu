function mesh = scale_mesh(mesh, scale, base)
%SCALE_MESH Scale mesh (NiHu / meshing)
%   MESH = SCALE_MESH(MESH, SCALE) scales the mesh MESH using the
%   scale factor(s) given by the scalar or vector SCALE
%
%   MESH = SCALE_MESH(MESH, SCALE, BASE) scales with respect to the vector
%   BASE
%
% See also: TRANSLATE_MESH, ROTATE_MESH, REVOLVE_MESH, EXTRUDE_MESH,
% REPEAT_MESH, REFLECT_MESH

%   Copyright 2008-2012 P. Fiala, P. Rucz
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modifed: 2012.12.14.

% process arguments
switch length(scale)
    case 1
        T = scale;
    case 2
        T = diag([scale(:)' 1]);
    case 3
        T = diag(scale);
    otherwise
        error('NiHu:scale_mesh:argFormat', ...
            'Invalid scaling size (%dx%d)', size(scale,1), size(scale,2));
end

% scale
if nargin > 2
    mesh = translate_mesh(mesh, -base);
end
mesh.Nodes(:,2:4) = mesh.Nodes(:,2:4) * T;
if nargin > 2
    mesh = translate_mesh(mesh, base);
end

end
