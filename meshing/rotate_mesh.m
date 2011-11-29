function mesh = rotate_mesh(mesh, varargin)
%ROTATE_MESH Rotate mesh around a given vector
%   MESH = ROTATE_MESH(MESH, DIR, DEG) or
%   MESH = ROTATE_MESH(MESH, BASE, DIR, DEG) rotates the mesh MESH
%   along the central line specified by BASE and DIR by and amount of
%   DEG (expressed in radians).
%
% See also: TRANSLATE_MESH, SCALE_MESH, REVOLVE_MESH, EXTRUDE_MESH,
% REPEAT_MESH, REFLECT_MESH

%   Copyright 2008-2010 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modifed: 02.12.2009

%% Argument check
error(nargchk(3, 4, nargin, 'struct'));
switch nargin
    case 3
        base = [0 0 0];
        dir = varargin{1};
        phi = varargin{2};
    case 4
        base = varargin{1};
        dir = varargin{2};
        phi = varargin{3};
end

coord = mesh.Nodes(:,2:4);
nNod = size(coord,1);
dir = dir/norm(dir); % unit normal vector
d = repmat(dir,nNod,1);
b = repmat(base,nNod,1);
x = coord - b;
a = repmat(dot(x,d,2),1,3).*d;
w1 = x - a;
w2 = cross(d,w1,2);
mesh.Nodes(:,2:4) = (a + cos(phi)*w1+sin(phi)*w2)+b;
end
