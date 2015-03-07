function [nodes, elements] = extract_core_mesh(model)
%EXTRACT_CORE_MESH Extract C++ core bem mesh from a NiHu mesh structure

%% Search for QUAD or TRIA elements and allocate space for bem model
ind12 = find(model.Elements(:,2) == ShapeSet.LinearLine.Id); % line elements
ind3 = find(model.Elements(:,2) == ShapeSet.LinearTria.Id); % linear tri elements
ind4 = find(model.Elements(:,2) == ShapeSet.LinearQuad.Id); % linear quad elements
ind23 = find(model.Elements(:,2) == ShapeSet.QuadraticTria.Id); % quadratic tria elements
ind24 = find(model.Elements(:,2) == ShapeSet.QuadraticQuad.Id); % quadratic quad elements
elem = drop_IDs(model);
elements = zeros(size(model.Elements,1),9+1);

%% build new element matrix
if ~isempty(ind12)
    elements(ind12,1) = 21202;
    elements(ind12,2:3) = elem(ind12,4+(1:2))-1; % C-style indexing
end
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
%     elements(ind24,1) = 32409;
    elements(ind24,1) = 32408;
    elements(ind24,1+(1:8)) = elem(ind24,4+(1:8))-1; % C-style indexing
end

%% build new node matrix
nodes = model.Nodes(:,2:end);
end
