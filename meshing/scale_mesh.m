function mesh = scale_mesh(mesh, scale)
%SCALE_MESH Scale mesh
%   MESH = SCALE_MESH(MESH, SCALE) scales the mesh MESH using the
%   scale factor(s) given by the scalar(vector) SCALE
%
% See also: TRANSLATE_MESH, ROTATE_MESH, REVOLVE_MESH, EXTRUDE_MESH,
% REPEAT_MESH, REFLECT_MESH

%   Copyright 2008-2010 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modifed: 02.12.2009

%% Argument check
error(nargchk(2, 2, nargin, 'struct'));

%%
switch length(scale)
    case 1
        mesh.Nodes(:,2:4) = mesh.Nodes(:,2:4) * scale;
    otherwise
        ind = 1+(1:length(scale));
        mesh.Nodes(:,ind) = mesh.Nodes(:,ind) .* ...
            repmat(scale(:).',size(mesh.Nodes,1),1);
end
end
