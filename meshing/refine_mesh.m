function newmesh = refine_mesh(mesh)
%REFINE_MESH refine mesh by dividing its elements into subelements
% NEWMESH = REFINE_MESH(MESH) refines the mesh by dividing its elements
% into sublements. Division is performed so that the element length is
% halved.
%
% Example:
%  mesh = create_sphere_boundary(1, 5);
%  m = refine_mesh(mesh);
%  figure; plot_mesh(mesh);
%  figure; plot_mesh(m);

elements = drop_IDs(mesh);

elem_data = {
    %ID, num_new_elems, num_new_nod, refine_func
    12, 2, 3, @refine_line
    23, 4, 6, @refine_tria
    24, 4, 9, @refine_quad
    };

newmesh.Elements = [];
newmesh.Nodes = [];

for iType = 1 : size(elem_data,1)
    id = elem_data{iType,1};
    nelems = elem_data{iType,2};    % num of new elements
    nnodes = elem_data{iType,3};    % new of new nodes
    refine_fun = elem_data{iType,4};    % refine function
    norignodes = mod(id,10);    % num of original nodes
    
    elems = elements(elements(:,2) == id, :); % select appropriate elements
    nel = size(elems,1);
    newelemnodes = zeros(nnodes*nel,3);
    newelemelements = zeros(nelems*nel,norignodes);
    for q = 1 : nel
        [newnod, newelem] = refine_fun(mesh.Nodes(elems(q,4+(1:norignodes)),2:end));
        newelemelements((q-1)*nelems+(1:nelems),1:norignodes) = newelem + (q-1)*nnodes;
        newelemnodes((q-1)*nnodes+(1:nnodes),1:3) = newnod;
    end
    
    newmesh.Elements(end+(1:size(newelemelements,1)),2:(4+norignodes)) =...
        [repmat([id 1 1], size(newelemelements,1), 1) newelemelements+size(newmesh.Nodes,1)];
    newmesh.Nodes(end+(1:size(newelemnodes,1)),2:4) = newelemnodes;
end

newmesh.Nodes(:,1) = 1 : size(newmesh.Nodes,1);
newmesh.Elements(:,1) = 1 : size(newmesh.Elements,1);

newmesh.Properties = mesh.Properties;
newmesh.Materials = mesh.Materials;

newmesh = merge_coincident_nodes(newmesh);

end

function [nodes, elements] = refine_line(nodes)
%REFINE_LINE refine a single line element into subelements
%   [NODES, ELEMENTS] = REFINE_LINE(NODES)

nodes = [
    nodes
    mean(nodes,1)
    ];
elements = [
    1 3
    3 2
    ];
end

function [nodes, elements] = refine_tria(nodes)
%REFINE_TRIA refine a single tria element into subelements
%   [NODES, ELEMENTS] = REFINE_TRIA(NODES)

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
%REFINE_QUAD refine a single quad element into subelements
%   [NODES, ELEMENTS] = REFINE_QUAD(NODES)

nodes = [
    nodes
    (nodes(1:4,:) + nodes([2 3 4 1], :)) /2
    mean(nodes,1)
    ];
elements = [
    1 5 9 8
    2 6 9 5
    3 7 9 6
    4 8 9 7
    ];
end
