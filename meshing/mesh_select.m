function [nodes elements] = mesh_select(mesh, expression, indmode, selmode)
%MESH_SELECT Element and node selection
%   [NODES, ELEMENTS] = MESH_SELECT(MESH, COND, INDMODE, SELMODE) selects
%   nodes and elements satisfying a given condition.
%   MESH : NiHU mesh structure
%   COND : list of selected node IDs or
%          any valid Matlab expression string containing the variables
%          'x', 'y', 'z', 'r', 'phi', 'R'. These variables refer to columns
%          of the matrix MESH.Nodes as follows:
%          'x'   = MESH.Nodes(:,2)
%          'y'   = MESH.Nodes(:,3)
%          'z'   = MESH.Nodes(:,4)
%          'r'   = sqrt(x^2 + y^2), 2D radius
%          'phi' = atan2(y, x), 2D angle
%          'theta' = atan2(z, r), 3D inclination angle
%          'R'   = sqrt(x^2 + y^2 + z^2), 3D radius
%   INDMODE : 'ID'  -> Node and Element IDs are returned (default)
%             'ind' -> Node and Element indices are returned. These indices
%                      index the matrices MESH.Nodes and MESH.Elements
%   SELMODE : 'all' -> An element is selected if all of its nodes satisfy
%                      the condition (default).
%             'any' -> An element is selected if any of its nodes satisfy
%                      the condition.
%             'none'-> An element is selected if none of its nodes satisfy
%                      the condition.
%             'nall'-> An element is selected if not all of its nodes
%                      satisfy the condition.
%
%   Example:
%      sphere = create_sphere_boundary(1,10);
%      [nodind, elind] = mesh_select(sphere, 'x+y < 1e-3', 'ind', 'all');
%      sphere.Elements = sphere.Elements(elind,:);
%      sphere = drop_unused_nodes(sphere);
%      plot_mesh(sphere);
%
% See also: MESH_SECTION, DROP_UNUSED_NODES

%   Copyright 2008-2010 P. Fiala, P. Rucz
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last updated: 06.09.2010.

%% Parameter check
error(nargchk(2, 4, nargin, 'struct'));
if nargin < 4
    selmode = 'all';
end
if nargin < 3
    indmode = 'ID';
end

%% Defining the expression variables
x = mesh.Nodes(:,2);
y = mesh.Nodes(:,3);
z = mesh.Nodes(:,4);
r = sqrt(x.^2 + y.^2); 
theta = atan2(z,r); %#ok<NASGU>
R = sqrt(x.^2 + y.^2 + z.^2); %#ok<NASGU>
phi = atan2(y,x); %#ok<NASGU>

%% Expression evaluation
Elements = drop_IDs(mesh);
nElements = size(Elements,1);
% find node indices satisfying the condition
if ischar(expression)
    nodeind = find(eval(expression));
else
    [trash, nodeind] = ismember(expression, mesh.Nodes(:,1)); %#ok<ASGLU>
end
switch lower(selmode)
    case 'any'
        elemind = find(any(ismember(Elements(:,5:end), nodeind),2));
    case 'none'
        elemind = find(any(ismember(Elements(:,5:end), nodeind),2));
        elemind = setdiff(1:nElements,elemind);
    case 'all'
        elemind = find(sum(ismember(Elements(:,5:end), nodeind),2) ==...
            sum(Elements(:,5:end) ~=0, 2));
    case 'nall'
        elemind = find(sum(ismember(Elements(:,5:end), nodeind),2) ==...
            sum(Elements(:,5:end) ~=0, 2));
        elemind = setdiff(1:nElements,elemind);
    otherwise
        error('Invalid selection mode selector: %s', selmode);
end
% IDs
nodeID = mesh.Nodes(nodeind, 1);
elemID = mesh.Elements(elemind, 1);

%% output selection
switch lower(indmode)
    case 'id'
        nodes = nodeID;
        elements = elemID;
    case 'ind'
        nodes = nodeind;
        elements = elemind;
    otherwise
        error('Invalid index mode selector: %s', indmode);
end
end
