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
%   MESH : AcouFEM mesh structure
%   DIR   : 1x3 vector corresponding to one segment of the extrudation
%   NREP  : numer of repetitions
%
% See also: TRANSLATE_MESH, SCALE_MESH, ROTATE_MESH, REVOLVE_MESH,
% REPEAT_MESH, REFLECT_MESH

%   Copyright 2008-2012 P. Fiala, P. Rucz
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modifed: 2012.12.19.

dir = dir(:).'; % ensure that dir is a row vector

% Create new nodes
nNod = size(mesh.Nodes,1);
coord = zeros((nRep+1)*nNod,3);
for iRep = 0 : nRep
    coord(iRep*nNod+(1:nNod),:) = bsxfun(@plus, mesh.Nodes(:,2:4), iRep*dir);
end

% Create new elements
Elements = drop_IDs(mesh);
% zero padding
Elements = [Elements zeros(size(Elements,1), 8-size(Elements,2))];

% Get cell normals
[~, normal] = centnorm(mesh);

% Select which elements to flip to conserve outward normals
reverse = (normal * dir.') < 0 & Elements(:,2) ~= ShapeSet.LinearLine.Id;

quads = Elements(:,2) == ShapeSet.LinearQuad.Id;
quad = Elements(quads,:);
nq = size(quad,1);
trias = Elements(:,2) == ShapeSet.LinearTria.Id;
tria = Elements(trias,:);
nt = size(tria,1);
lines = Elements(:,2) == ShapeSet.LinearLine.Id;
line = Elements(lines,:);
nl = size(line,1);
nE = nq+nt+nl;

flipquad = reverse(quads);
fliptria = reverse(trias);
flipline = reverse(lines);

quad(flipquad,5:8) = fliplr(quad(flipquad,5:8));
tria(fliptria,5:7) = fliplr(tria(fliptria,5:7));
line(flipline,5:6) = fliplr(line(flipline,5:6));

Elem = zeros(nRep*nE,12);

for iRep = 1 : nRep
    Elem((iRep-1)*nE+(1:nq),1:12) = [repmat([0 ShapeSet.LinearHexa.Id],nq,1) quad(:,3:4) quad(:,5:8)+(iRep-1)*nNod quad(:,5:8)+iRep*nNod];
    Elem((iRep-1)*nE+nq+(1:nt),1:10) = [repmat([0 ShapeSet.LinearPenta.Id],nt,1) tria(:,3:4) tria(:,5:7)+(iRep-1)*nNod tria(:,5:7)+iRep*nNod];
    Elem((iRep-1)*nE+nq+nt+(1:nl),1:8) = [repmat([0 ShapeSet.LinearQuad.Id],nl,1) line(:,3:4) line(:,[5 6])+(iRep-1)*nNod line(:,[6 5])+iRep*nNod];
end

% Assemble new mesh structure
mesh2.Nodes(:,2:4) = coord;
mesh2.Nodes(:,1) = 1:size(mesh2.Nodes,1);
mesh2.Elements = Elem;
mesh2.Elements(:,1) = 1:size(mesh2.Elements,1);
mesh2.Properties = mesh.Properties;
mesh2.Materials = mesh.Materials;
end
