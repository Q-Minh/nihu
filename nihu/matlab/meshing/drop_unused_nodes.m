function [mesh, nodind, dropind] = drop_unused_nodes(mesh)
%DROP_UNUSED_NODES  Exclude unused nodes from a mesh
%   [MESH, NODIND] = DROP_UNUSED_NODES(MESH) returns a new MESH,
%   where the matrix MESH.Nodes does not contain the nodes that are not
%   referred to by the matrix MESH.Elements. The optional output parameter
%   NODIND contains the indices of the used nodes in the original MESH
%   structure, DROPIND contains the indices of the nodes that have been
%   dropped.
%
% See also: MERGE_COINCIDENT_NODES, DROP_IDS

%   Copyright 2008-2010 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last updated: 2012.12.19.

nodind = find(ismember(mesh.Nodes(:,1), mesh.Elements(:,5:end)));
if nargout > 2
    dropind = setdiff(1:size(mesh.Nodes,1), nodind);
end
mesh.Nodes = mesh.Nodes(nodind,:);
end
