function I = cell2node_interp(model)
%CELL2NODE_INTERP Summary of this function goes here
%   Detailed explanation goes here

% ONLY WORKS FOR TRIA MESH!!

elem = drop_IDs(model);
nElem = size(elem,1);
nNodes = size(model.Nodes,1);

N = 4;

% Create node vector for each elem
allNode = elem(:,4+(1:N)).';
allNode = allNode(:);

% Create elem ID vector (N x 1 IDs under each other)
allElem = repmat(elem(:,1),1,N).';
allElem = allElem(:);

allValues = ones(size(allNode));

[allNode, ind] = sort(allNode, 'ascend');
allElem = allElem(ind);

cStart = 1;
cPos = 1;
cNode = allNode(1);

while cPos < length(allNode);
    if allNode(cPos+1) ~= cNode
        allValues(cStart:cPos) = allValues(cStart:cPos)/(cPos-cStart+1);
        cStart = cPos+1;
        cNode = allNode(cPos+1);
    end
    cPos = cPos +1;
end

I = sparse(allNode, allElem, allValues, nNodes, nElem);

end

