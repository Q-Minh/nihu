function mesh2 = revolve_mesh(mesh, ang, nRep, dir, base)
%REVOLVE_MESH Revolve 1D and 2D mesh around a given vector (NiHu / meshing)
%   MESH = REVOLVE_MESH(MESH, ANG, NREP, DIR) or
%   MESH = REVOLVE_MESH(MESH, ANG, NREP, DIR, BASE) revolves the 1D or 2D
%   NiHu mesh mesh around the central line given by the vectors BASE and
%   DIR, and creates a 2D or 3D mesh.
%   The initial 1D mesh contains LINE elements, and the resulting mesh
%   will contain QUAD elements.
%   The TRIA and QUAD elements of the initial 2D mesh will be transformed
%   into PENTA and HEXA elements, respectively.
%   The geometry is repeated NPHI times, the angular increment of the
%   repetitions is DPHI.
%
% See also: TRANSLATE_MESH, SCALE_MESH, ROTATE_MESH, EXTRUDE_MESH,
% REPEAT_MESH, REFLECT_MESH

%   Copyright 2008-2012 P. Fiala, P. Rucz
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modifed: 2015.03.16.

% translate mesh to base if needed
if nargin > 4
    mesh = translate_mesh(mesh, -base);
end

% compute matrix of a single rotation
T = rotation_matrix(ang, dir);

% Create new nodes
coord = mesh.Nodes(:,2:4);
nNod = size(coord,1);

coords = zeros((nRep+1)*nNod,3);
for iPhi = 0 : nRep
    coords(iPhi*nNod+(1:nNod),:) = coord;
    coord = coord * T;
end

% Select which elements to flip to conserve outward normals
[~, normal] = centnorm(mesh);
lsetid = mesh.Elements(:,2);
reverse = false(size(mesh.Elements,1),1);

Elem = extrude_elements(mesh, nNod, nRep, reverse);

% Assemble new mesh structure
mesh2.Nodes(:,2:4) = coords;
mesh2.Nodes(:,1) = 1:size(mesh2.Nodes,1);
mesh2.Elements = Elem;
mesh2.Elements(:,1) = 1:size(mesh2.Elements,1);
mesh2.Properties = mesh.Properties;
mesh2.Materials = mesh.Materials;

% translate back if needed
if nargin > 4
    mesh2 = translate_mesh(mesh2, base);
end

end
