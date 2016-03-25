function mesh2 = extrude_mesh(mesh, dir, nRep)
%EXTRUDE_MESH Extrude 1D and 2D mesh along a given direction (NiHu / meshing)
%   MESH = EXTRUDE_MESH(MESH, DIR, NREP) extrudes the 1D or 2D
%   FE mesh MESH along the direction DIR and creates a new FE mesh.
%   The LINE elements of an initial 1D mesh are extruded into QUAD elements.
%   The TRIA and QUAD elements of the initial 2D mesh are extruded into
%   PENTA and HEXA elements, respectively.
%   The geometry is repeated NREP times.
%
%   Parameters:
%   MESH : NiHu mesh structure
%   DIR   : 1x3 vector corresponding to one segment of the extrudation
%   NREP  : numer of repetitions
%
% See also: TRANSLATE_MESH, SCALE_MESH, ROTATE_MESH, REVOLVE_MESH,
% REPEAT_MESH, REFLECT_MESH

%   Copyright 2008-2012 P. Fiala, P. Rucz
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modifed: 2015.03.16.

dir = dir(:).'; % ensure that dir is a row vector

% Create new nodes
nNod = size(mesh.Nodes,1);
coord = zeros((nRep+1)*nNod,3);
for iRep = 0 : nRep
    coord(iRep*nNod+(1:nNod),:) = bsxfun(@plus, mesh.Nodes(:,2:4), iRep*dir);
end

% Select which elements to flip to conserve outward normals
[~, normal] = centnorm(mesh);
lsetid = mesh.Elements(:,2);
reverse = lsetid ~= ShapeSet.LinearLine.Id & (normal * dir.') < 0 | ...
    lsetid == ShapeSet.LinearLine.Id & (normal * dir.') > 0;
mesh.Elements(reverse,:) = flip_elements(mesh.Elements(reverse,:));

Elem = drop_IDs(mesh);
Elem = extrude_elements(Elem, nNod, nRep);

% Assemble new mesh structure
mesh2.Nodes(:,2:4) = coord;
mesh2.Nodes(:,1) = 1:size(mesh2.Nodes,1);
mesh2.Elements = Elem;
mesh2.Elements(:,1) = 1:size(mesh2.Elements,1);
mesh2.Properties = mesh.Properties;
mesh2.Materials = mesh.Materials;

end
