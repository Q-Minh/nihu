function v = node_volume(model)
% Calculates corresponding volume for nodes

elem = drop_IDs(model);

elem_v = mesh_volume(model);

v = zeros(size(model.Nodes,1),1);

for iElem = 1:size(elem,1)
    type = elem(iElem,2);
    nodenum = mod(type,10);
    nodinds = elem(iElem,4+(1:nodenum));
    v(nodinds) = v(nodinds) + 1/nodenum*elem_v(iElem);
end

end