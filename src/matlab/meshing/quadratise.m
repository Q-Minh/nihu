function newmesh = quadratise(mesh, options)
%QUADRATISE replace linear elements with quadratic elements

if nargin < 2
    options = struct([]);
end

elements = drop_IDs(mesh);

elem_data = {
    %ID, num_new_elems, num_new_nod, refine_func
    ShapeSet.LinearLine.Id, 2, 3, @quadratise_line
    ShapeSet.LinearTria.Id, 4, 6, @quadratise_tria
    ShapeSet.LinearQuad.Id, 4, 9, @quadratise_quad
};

newmesh.Properties = mesh.Properties;
newmesh.Materials = mesh.Materials;
newmesh.Elements = [];
newmesh.Nodes = [];

for iType = 1 : size(elem_data,1)
    id = elem_data{iType,1};
    nelems = elem_data{iType,2};    % num of new elements
    nnodes = elem_data{iType,3};    % new of new nodes
    refine_fun = elem_data{iType,4};    % refine function
    norignodes = size(ShapeSet.fromId(id).Nodes,1);    % num of original nodes
    
    elems = elements(elements(:,2) == id, :); % select appropriate elements
    if isempty(elems)
        continue;
    end
    nel = size(elems,1);
    newelemnodes = zeros(nnodes*nel,3);
    newelemelements = zeros(nel,norignodes);
    for q = 1 : nel
        [newnod, newelem] = refine_fun(mesh.Nodes(elems(q,4+(1:norignodes)),2:end), options);
        newelemelements(q,1:length(newelem)) = newelem + (q-1)*nnodes;
        newelemnodes((q-1)*nnodes+(1:nnodes),1:3) = newnod;
    end
    
    newmesh.Elements(end+(1:size(newelemelements,1)),2:(4+size(newelemelements,2))) =...
        [repmat([quadratise_id(id, options) 1 1], size(newelemelements,1), 1) newelemelements+size(newmesh.Nodes,1)];
    newmesh.Nodes(end+(1:size(newelemnodes,1)),2:4) = newelemnodes;
end

newmesh.Nodes(:,1) = 1 : size(newmesh.Nodes,1);
newmesh.Elements(:,1) = 1 : size(newmesh.Elements,1);

warning('TODO: merge coincident nodes does not work well with quadratic elements');
newmesh = merge_coincident_nodes(newmesh, 1e-3);

end


function newid = quadratise_id(id, options)
    switch id
        case ShapeSet.LinearLine.Id 
            newid = ShapeSet.QuadraticLine.Id;
        case ShapeSet.LinearTria.Id
            newid = ShapeSet.QuadraticTria.Id;
        case ShapeSet.LinearQuad.Id
            if (isfield(options, 'quad_with_mid') && options.quad_with_mid)
                newid = ShapeSet.QuadraticQuadMid.Id;
            else
                newid = ShapeSet.QuadraticQuad.Id;
            end
    end
%newid = 10 * id + 2;
%newid = 10000 + id;
end

function [nodes, element] = quadratise_line(nodes, options)

nodes = [
    nodes
    mean(nodes,1)
    ];
element = [1 3 2];
end

function [nodes, element] = quadratise_tria(nodes, options)

nodes = [
    nodes
    (nodes(1:3,:) + nodes([2 3 1], :)) /2
    ];
element = [1 4 2 5 3 6];
end

function [nodes, element] = quadratise_quad(nodes, options)

if isfield(options, 'quad_with_mid') && options.quad_with_mid
    nodes = [
        nodes
        (nodes(1:4,:) + nodes([2 3 4 1], :)) /2
        mean(nodes,1)
        ];
    element = [1 5 2 6 3 7 4 8 9];

else
    nodes = [
        nodes
        (nodes(1:4,:) + nodes([2 3 4 1], :)) /2
    ];
    element = [1 5 2 6 3 7 4 8];
end
    
end
