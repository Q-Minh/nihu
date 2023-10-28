function d = mesh_edge_size(model)
%MESH_EDGE_SIZE Estimate edge size of a mesh (NiHu / meshing)
%   D = MESH_EDGE_SIZE(MESH) returns the estimated edge size for all
%   elements of the MESH.

%   Copyright 2009-2015 P. Fiala, P. Rucz
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last updated: 2015.04.01. FP. edge rules introduced

% preproc
elem = drop_IDs(model);
d = zeros(size(elem,1),1);

edgeRules = {
    ShapeSet.LinearLine.Id, [1 2]
    ShapeSet.LinearTria.Id, [1 2; 2 3; 3 1]
    ShapeSet.LinearQuad.Id, [1 2; 2 3; 3 4; 4 1]
    ShapeSet.LinearTetra.Id, [1 2; 2 3; 3 1; 1 4; 2 4; 3 4]
    ShapeSet.LinearPenta.Id, [1 2; 2 3; 3 1; 1 4; 2 5; 3 6; 4 5; 5 6; 6 4]
    ShapeSet.LinearHexa.Id, [1 2; 2 3; 3 4; 4 1; 1 5; 2 6; 3 7; 4 8; 5 6; 6 7; 7 8; 8 5]
    };

srcids = cell2mat(edgeRules(:,1));

uniqueIds = unique(elem(:,2));
for i = 1 : length(uniqueIds)
    lsetid = uniqueIds(i);
    sel = elem(:,2) == lsetid;
    q = lsetid == srcids;
    if ~any(q)
        error('NiHu:mesh_edge_size:TODO', 'lset id %d not supported yet', lsetid);
    end
    edges = edgeRules{q,2};
    
    for k = 1 : size(edges,1)
        evec = model.Nodes(elem(sel,4+edges(k,2)),2:end)-...
            model.Nodes(elem(sel,4+edges(k,1)),2:end);
        e = sqrt(dot(evec,evec,2));
        d(sel) = max(d(sel), e);
    end
end
end
