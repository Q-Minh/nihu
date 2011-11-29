function meshes = split_independent_meshes(mesh)
%SPLIT_INDEPENDENT_MESHES Split mesh to independent submeshes
%   MESHES = SPLIT_INDEPENDENT_MESHES(MESH)  splits the mesh MESH into
%   independent parts. Independent parts don't share nodes. MESHES is an
%   array of independent mesh structures.
%
% See also: JOIN_MESHES

%   Copyright 2008-2010 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modifed: 02.12.2009

%% Argument check
error(nargchk(1, 1, nargin, 'struct'));

%% get rid of IDs
elem = drop_IDs(mesh);
elem = elem(:,5:end);

%% start at node 1 and find connected nodes level-by-level
node0 = 1; % index of starting node
while true
    % find neighbors
    nodes = elem(any(ismember(elem, node0),2),:);
    nodes = unique(nodes(nodes~=0));
    % if no new neighbors have been found, stop
    if length(nodes) == length(node0)
        break;
    else
        node0 = nodes;
    end
end

%% The first mesh consists of the newly found part
meshes.Nodes = mesh.Nodes(nodes,:);
meshes.Elements = mesh.Elements(any(ismember(elem, nodes), 2), :);
meshes.Properties = mesh.Properties;
meshes.Materials = mesh.Materials;

%% The second mesh is the rest, and is recursively further divided
rest = setdiff(1:size(mesh.Nodes,1), nodes);
if ~isempty(rest)
    mesh2.Nodes = mesh.Nodes(rest,:);
    mesh2.Elements = mesh.Elements(any(ismember(elem, rest), 2), :);
    mesh2.Properties = mesh.Properties;
    mesh2.Materials = mesh.Materials;
    meshes = [meshes split_independent_meshes(mesh2)];
end

end
