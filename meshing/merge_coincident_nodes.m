function [mesh, nodind, dropind] = merge_coincident_nodes(mesh, tol)
%MERGE_COINCIDENT_NODES Merge coincident nodes in a mesh
%   [MESH2, NODIND, DROPIND] = MERGE_COINCIDENT_NODES(MESH2, TOL) searches
%   and merges the coincident nodes in the MESH. Nodes are coincident if
%   their distance is less than TOL. The default value for TOL is 1e-3.
%   The output parameter NODIND contains the indices of the nodes that are
%   contained in the output mesh structure MESH2 so that MESH2.Nodes(:,1) =
%   MESH.Nodes(nodind,1). The output parameter DROPIND contains the indices
%   of the nodes that have been dropped.
%
%   Example:
%    cross = create_slab(1,10);
%    cross = translate_mesh(cross, [1 0 0]);
%    vol = revolve_mesh(cross, [0 0 0], [0 1 0], 2*pi/40, 40);
%    [vol2, nodind, dropind] = merge_coincident_nodes(vol, 1e-2);
%    fprintf(1, '%d nodes merged\n', length(dropind));
%    plot_fem_model(vol2);
%
% See also: DROP_UNUSED_NODES, DROP_IDS

%   Copyright 2008-2010 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last updated: 02.12.2009.

%% Parameter check
error(nargchk(1, 2, nargin, 'struct'));
if nargin < 2
    tol = 1e-3;
end

%% get rid of original IDs
elements = drop_IDs(mesh);

%% Coincident nodes
[trash, nodind, perm] = unique(round(mesh.Nodes(:,2:4)/tol), 'rows');
if nargout > 2
    dropind = setdiff(1:size(mesh.Nodes,1), nodind);
end
mesh.Nodes = mesh.Nodes(nodind,:);

%% Renumber the node IDs
p = [0; mesh.Nodes(perm,1)];
mesh.Elements(:,5:end) = p(elements(:,5:end)+1);
end
