function [nodes, elements] = extract_core_mesh(model)
%EXTRACT_CORE_MESH Extract C++ core bem mesh from a NiHu mesh structure

elem = drop_IDs(model);
elements = zeros(size(model.Elements,1),9+1);

data = {
    ShapeSet.LinearLine, 21202
    ShapeSet.LinearTria, 32303
    ShapeSet.LinearQuad, 32404
    ShapeSet.QuadraticTria, 32306
    ShapeSet.QuadraticQuad, 32408
    ShapeSet.QuadraticQuadMid, 32409
    };

for i = 1 : size(data,1)
    lset = data{i,1};
    core_id = data{i,2};
    ind = find(model.Elements(:,2) == lset.Id);
    if ~isempty(ind)
        nNode = size(lset.Nodes,1);
        elements(ind,1) = core_id;
        elements(ind,1+(1:nNode)) = elem(ind,4+(1:nNode))-1;
    end
end

nodes = model.Nodes(:,2:end);
end
