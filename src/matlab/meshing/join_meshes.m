function mesh = join_meshes(mesh1, mesh2, varargin)
%JOIN_MESHES  Create a single mesh from several submeshes
%   MESH = JOIN_MESHES(MES1, MESH2, ...) Combines the meshes MESH1,
%   MESH2, ... into a single mesh structure. The resulting structure
%   contains all the nodes, properties and materials of all the meshes,
%   without merging the duplicated entries. The first lines of the
%   matrix MESH.Nodes correspond to MESH1.Nodes, the second block
%   corresponds to MESH2.Nodes, etc., without resorting the entries. The
%   same is true for the matrices MESH.Properties, Materials and Elements.
%
% See also: SPLIT_INDEPENDENT_MESHES, MERGE_COINCIDENT_NODES

%   Copyright 2008-2012 P. Fiala, P. Rucz
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modifed: 2012.12.19.

% Recursive mechanism
if nargin > 2
    mesh = join_meshes(join_meshes(mesh1, mesh2), varargin{:});
else
    % We'll have new IDs
    e1 = drop_IDs(mesh1);
    e2 = drop_IDs(mesh2);
    
    % Material and properties compatibility
    s1 = size(mesh1.Materials,2);
    s2 = size(mesh2.Materials,1);
    if (s1 > s2)
        mesh2.Materials(:,s1) = 0;
    elseif s2 > s1
        mesh1.Materials(:,s2) = 0;
    end
    mesh.Materials = union(mesh1.Materials(:,2:end),mesh2.Materials(:,2:end),'rows');
    [~, matInd1, matOld1] = intersect(mesh.Materials,mesh1.Materials(:,2:end),'rows');
    [~, matInd2, matOld2] = intersect(mesh.Materials,mesh2.Materials(:,2:end),'rows');
    mesh.Materials = [(1:size(mesh.Materials,1)).',mesh.Materials];

    s1 = size(mesh1.Properties,2);
    s2 = size(mesh2.Properties,1);
    if (s1 > s2)
        mesh2.Properties(:,s1) = 0;
    elseif s2 > s1
        mesh1.Properties(:,s2) = 0;
    end
    Prop = union(mesh1.Properties(:,2:end),mesh2.Properties(:,2:end),'rows');
    [~, propInd1, propOld1] = intersect(Prop,mesh1.Properties(:,2:end),'rows');
    [~, propInd2, propOld2] = intersect(Prop,mesh2.Properties(:,2:end),'rows');
    mesh.Properties = [(1:size(Prop,1)).',Prop];
    
    % Get number of nodes and elements
    nN1 = size(mesh1.Nodes,1);  % number of nodes in meshes
    nN2 = size(mesh2.Nodes,1);
    nE1 = size(e1,1);           % number of elements in meshes
    nE2 = size(e2,1);
    
    % Assign nodes
    mesh.Nodes = zeros(nN1+nN2,4);
    mesh.Nodes(1:nN1,:) = mesh1.Nodes;
    mesh.Nodes(nN1+1:end,:) = mesh2.Nodes;
    mesh.Nodes(:,1) = 1:(nN1+nN2);
    
    % Assign elements
    mesh.Elements = zeros(nE1+nE2,max([size(e1,2), size(e2,2)]));
    mesh.Elements(1:nE1,1:size(e1,2)) = e1;
    elem = e2(:,5:end);
    elem(elem~=0) = elem(elem~=0) + nN1;
    mesh.Elements(nE1+(1:nE2),4+(1:size(elem,2))) = elem;
    
    % Reassign materials and properties
    mesh.Elements(1:nE1,1:4) = [e1(:,1:2),...
        mesh.Materials(matInd1(matOld1(e1(:,3))),1),...
        mesh.Properties(propInd1(propOld1(e1(:,4))),1)];
    mesh.Elements(nE1+(1:nE2),1:4) = [e2(:,1:2),...
        mesh.Materials(matInd2(matOld2(e2(:,3))),1),...
        mesh.Properties(propInd2(propOld2(e2(:,4))),1)];
    
    % Renumber elements
    mesh.Elements(:,1) = 1:(nE1+nE2);
end
