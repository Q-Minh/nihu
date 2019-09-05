function [nodes, elements] = extract_core_mesh(model, type)
%EXTRACT_CORE_MESH Extract C++ core bem mesh from a NiHu mesh structure

if nargin < 2
    type = 'surface';
end

elem = drop_IDs(model);
elements = zeros(size(model.Elements,1),9+1);

switch type
    case 'surface'
        data = {
            ShapeSet.LinearLine, 21202
            ShapeSet.LinearTria, 32303
            ShapeSet.LinearQuad, 32404
            ShapeSet.QuadraticTria, 32306
            ShapeSet.QuadraticQuad, 32408
            ShapeSet.QuadraticQuadMid, 32409
        };
    case 'volume'
        data = {
            ShapeSet.LinearLine, 11202
            ShapeSet.LinearTria, 22303
            ShapeSet.LinearQuad, 22404
            ShapeSet.QuadraticTria, 22306
            ShapeSet.QuadraticQuad, 22408
            ShapeSet.QuadraticQuadMid, 22409
        };
end

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
