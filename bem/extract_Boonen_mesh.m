function [nodes, elements] = extract_Boonen_mesh(model)
%EXTRACT_BOONEN_MESH Extract Boonen13 bem mesh from a NiHu mesh structure

%% Search for QUAD or TRIA elements and allocate space for bem model
ind3 = find(model.Elements(:,2) == 23); % tri elements
ind4 = find(model.Elements(:,2) == 24); % quad elements
ind23 = find(model.Elements(:,2) == 232); % quad elements
ind24 = find(model.Elements(:,2) == 242); % quad elements
elem = drop_IDs(model);
elements = zeros(size(model.Elements,1),9+1);

%% build new element matrix
if ~isempty(ind3)
    elements(ind3,1) = 32303;
    elements(ind3,2:4) = elem(ind3,4+(1:3))-1; % C-style indexing
end
if ~isempty(ind4)
    elements(ind4,1) = 32404;
    elements(ind4,2:5) = elem(ind4,4+(1:4))-1; % C-style indexing
end
if ~isempty(ind23)
    elements(ind23,1) = 32306;
    elements(ind23,1+(1:6)) = elem(ind23,4+(1:6))-1; % C-style indexing
end
if ~isempty(ind24)
    elements(ind24,1) = 32409;
    elements(ind24,2:10) = elem(ind24,4+(1:9))-1; % C-style indexing
end

%% build new node matrix
nodes = model.Nodes(:,2:4);
end
