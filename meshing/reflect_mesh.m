function mesh = reflect_mesh(mesh, varargin)
%REFLECT_MESH Reflect mesh to a symmetry plane
%   MESH = REFLECT_MESH(MESH, NORMAL) or 
%   MESH = REFLECT_MESH(MESH, BASE, NORMAL) reflects the NiHu mesh MESH to
%   the symmetry plane defined by BASE and NORMAL.
%   BASE   : a point of the symmetry plane, default is [0 0 0]
%   NORMAL : normal of the symmetry plane
%
% See also: TRANSLATE_MESH, SCALE_MESH, ROTATE_MESH, EXTRUDE_MESH,
% REVOLVE_MESH, REPEAT_MESH

%   Copyright 2008-2010 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last updated: 02.12.2009.

%% Argument check
error(nargchk(2, 3, nargin, 'struct'));
switch nargin
    case 2
        base = [0 0 0];
        normal = varargin{1};
    case 3
        base = varargin{1};
        normal = varargin{2};
end

%%
nN = size(mesh.Nodes,1);
base = repmat(base(:).', nN, 1);
normal = normal(:).';
normal = repmat(normal / norm(normal), nN, 1);

% reflect nodes
nodes = mesh.Nodes(:,2:4);
nodes = nodes - (2*repmat(dot(nodes-base, normal, 2),1,3)).*normal;
mesh.Nodes(:,2:4) = nodes;

% flip elements
mesh = flip_elements(mesh);
end
