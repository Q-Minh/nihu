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

% Last modifed: 2012.12.19.

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

% Create new elements
Elements = drop_IDs(mesh);
% zero padding
Elements = [Elements zeros(size(Elements,1), 8-size(Elements,2))];
quad = Elements(Elements(:,2) == 24,:);
nq = size(quad,1);
tria = Elements(Elements(:,2) == 23,:);
nt = size(tria,1);
line = Elements(Elements(:,2) == 12,:);
nl = size(line,1);
nE = nq+nt+nl;

Elem = zeros(nRep*nE,12);

for iPhi = 1 : nRep
    Elem((iPhi-1)*nE+(1:nq),1:12) = [repmat([0 38],nq,1) quad(:,3:4) quad(:,5:8)+(iPhi-1)*nNod quad(:,5:8)+iPhi*nNod];
    Elem((iPhi-1)*nE+nq+(1:nt),1:10) = [repmat([0 36],nt,1) tria(:,3:4) tria(:,5:7)+(iPhi-1)*nNod tria(:,5:7)+iPhi*nNod];
    Elem((iPhi-1)*nE+nq+nt+(1:nl),1:8) = [repmat([0 24],nl,1) line(:,3:4) line(:,[5 6])+(iPhi-1)*nNod line(:,[6 5])+iPhi*nNod];
end

% Assemble new mesh
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
