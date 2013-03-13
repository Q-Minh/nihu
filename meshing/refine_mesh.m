function newmesh = refine_mesh(mesh)
%REFINE_MESH refine mesh by dividing its tria and quad elements into subelements

elements = drop_IDs(mesh);

quads = elements(elements(:,2) == 24, :);
nquad = size(quads,1);
newquadnodes = zeros(9*nquad,3);
newquadelements = zeros(4*nquad,4);
for q = 1 : nquad
    [newnod, newelem] = refine_quad(mesh.Nodes(quads(q,5:8),2:4));
    newquadelements((q-1)*4+(1:4),1:4) = newelem + (q-1)*9;
    newquadnodes((q-1)*9+(1:9),1:3) = newnod;
end

trias = elements(elements(:,2) == 23, :);
ntria = size(trias,1);
newtrianodes = zeros(9*ntria,3);
newtriaelements = zeros(4*ntria,3);
for q = 1 : ntria
    [newnod, newelem] = refine_tria(mesh.Nodes(trias(q,5:7),2:4));
    newtriaelements((q-1)*4+(1:4),1:3) = newelem + (q-1)*6;
    newtrianodes((q-1)*6+(1:6),1:3) = newnod;
end


newmesh.Nodes(:,2:4) = [
    newquadnodes;
    newtrianodes;
    ];
newmesh.Nodes(:,1) = 1 : size(newmesh.Nodes,1);

newmesh.Elements(1:4*nquad,2:8) = [repmat([24 1 1], 4*nquad, 1) newquadelements];
newmesh.Elements(end+(1:4*ntria),2:7) = [repmat([23 1 1], 4*ntria, 1) newtriaelements+9*nquad];
newmesh.Elements(:,1) = 1 : size(newmesh.Elements,1);

newmesh.Properties = mesh.Properties;
newmesh.Materials = mesh.Materials;

newmesh = merge_coincident_nodes(newmesh);

end

function [nodes, elements] = refine_tria(nodes)
nodes = [
    nodes
    (nodes(1:3,:) + nodes([2 3 1], :)) /2
    ];
elements = [
    1 4 6
    2 5 4
    3 6 5
    4 5 6
    ];
end

function [nodes, elements] = refine_quad(nodes)
nodes = [
    nodes
    (nodes(1:4,:) + nodes([2 3 4 1], :)) /2
    sum(nodes,1) / 4
    ];
elements = [
    1 5 9 8
    2 6 9 5
    3 7 9 6
    4 8 9 7
    ];
end
