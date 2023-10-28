function meshes = split_independent_meshes(mesh)
%SPLIT_INDEPENDENT_MESHES Split mesh to independent submeshes
%   MESHES = SPLIT_INDEPENDENT_MESHES(MESH)  splits the mesh MESH into
%   independent parts. Independent parts don't share nodes. MESHES is an
%   array of independent mesh structures.
%
% See also: JOIN_MESHES

%   Copyright 2008-2012 P. Fiala, P. Rucz
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modifed: 2012.12.19.

[S, E] = adjacency(mesh);
part = node_partitions(S);
nParts = max(part);
meshes = struct([]);
for i = 1 : nParts
    meshes(i).Nodes = mesh.Nodes(part == i,:);
    meshes(i).Elements = mesh.Elements(any(E(:,part==i),2),:);
    meshes(i).Materials = mesh.Materials;
    meshes(i).Properties = mesh.Properties;
end

end
