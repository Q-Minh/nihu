function [nodes, elements] = extract_bem_model(model)
%extract_bem_model  Extract bem model from an acoufem model structure
%  [nodes, elements, nodid, elid] = extract_bem_model(model) extracts the
%  acoubem model fields from the acoufem model structure model. The new
%  model contains only the QUAD and TRIA elements of the original mesh.
% Parameters:
%  model    : AcouFEM model structure
% Output:
%  nodes    : nNx3 matrix containing xyz coordinates of surface vertices
%  elements : nEx5 matrix containing element vertex indices. Each row has
%             the structure
%             [n v1 v2 v3 v4]
%             where n = 3 or 4 for TRIA and QUAD elements, and vi refer to
%             the vertex indices in NODES. For TRIA elements, v4 = 0.

% Peter Fiala
% 2009

%% Search for QUAD or TRIA elements and allocate space for bem model
ind3 = find(model.Elements(:,2) == 23); % tri elements
ind4 = find(model.Elements(:,2) == 24); % quad elements
elem = drop_IDs(model);
elements = zeros(size(model.Elements,1),5);

%% build new element matrix
if ~isempty(ind3)
    elements(ind3,1) = 3;
    elements(ind3,2:4) = elem(ind3,5:7)-1;
end
if ~isempty(ind4)
    elements(ind4,1) = 4;
    elements(ind4,2:5) = elem(ind4,5:8)-1;
end

%% build new node matrix
nodes = model.Nodes(:,2:4);
end