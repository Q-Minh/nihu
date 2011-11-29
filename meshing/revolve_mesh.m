function mesh2 = revolve_mesh(mesh, varargin)
%REVOLVE_MESH Revolve 1D and 2D mesh around a given vector
%   MESH = REVOLVE_MESH(MESH, DIR, DPHI, NPHI) or
%   MESH = REVOLVE_MESH(MESH, BASE, DIR, DPHI, NPHI) revolves the 1D or 2D
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

%   Copyright 2008-2010 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modifed: 02.12.2009

%% Argument check
error(nargchk(4, 5, nargin, 'struct'));
switch nargin
    case 4
        base = [0 0 0];
        dir = varargin{1};
        dphi = varargin{2};
        nPhi = varargin{3};
    case 5
        base = varargin{1};
        dir = varargin{2};
        dphi = varargin{3};
        nPhi = varargin{4};
end

%% Create new nodes
coord = mesh.Nodes(:,2:4);
nNod = size(coord,1);
dir = dir./sqrt(dot(dir,dir)); % unit normal vector
d = repmat(dir,nNod,1);
b = repmat(base,nNod,1);
x = coord - b;
a = repmat(dot(x,d,2),1,3).*d;
w1 = x - a;
w2 = cross(d,w1,2);
coord = zeros(nPhi*nNod,3);

coord(1:nNod,:) = x+b;
for iPhi = 1 : nPhi
    phi = iPhi * dphi;
    coord(iPhi*nNod+(1:nNod),:) = (a + cos(phi)*w1+sin(phi)*w2)+b;
end

%% Create new elements
Elements = drop_IDs(mesh);
quad = Elements(Elements(:,2) == 24,:);
nq = size(quad,1);
tria = Elements(Elements(:,2) == 23,:);
nt = size(tria,1);
line = Elements(Elements(:,2) == 12,:);
nl = size(line,1);
nE = nq+nt+nl;

Elem = zeros(nPhi*nE,12);

for iPhi = 1 : nPhi
    if nq > 0
        Elem((iPhi-1)*nE+(1:nq),1:12) = [repmat([0 38],nq,1) quad(:,3:4) quad(:,5:8)+(iPhi-1)*nNod quad(:,5:8)+iPhi*nNod];
    end
    if nt > 0
        Elem((iPhi-1)*nE+nq+(1:nt),1:10) = [repmat([0 36],nt,1) tria(:,3:4) tria(:,5:7)+(iPhi-1)*nNod tria(:,5:7)+iPhi*nNod];
    end
    if nl > 0
        Elem((iPhi-1)*nE+nq+nt+(1:nl),1:8) = [repmat([0 24],nl,1) line(:,3:4) line(:,[5 6])+(iPhi-1)*nNod line(:,[6 5])+iPhi*nNod];
    end
end

%% Assemble new mesh
mesh2.Nodes(:,2:4) = coord;
mesh2.Nodes(:,1) = 1:size(mesh2.Nodes,1);
mesh2.Elements = Elem;
mesh2.Elements(:,1) = 1:size(mesh2.Elements,1);
mesh2.Properties = mesh.Properties;
mesh2.Materials = mesh.Materials;
end
