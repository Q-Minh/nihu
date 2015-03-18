function Elem = extrude_elements(Elements, nNod, nRep, opt)

if nargin < 4
    opt = 'finite';
end

if nargin < 3
    nRep = 1;
end

switch opt
    case 'finite'
        data = {
            ShapeSet.LinearLine, ShapeSet.LinearQuad, [1 2; 2 1]
            ShapeSet.LinearTria, ShapeSet.LinearPenta, [1 2 3; 1 2 3]
            ShapeSet.LinearQuad, ShapeSet.LinearHexa, [1 2 3 4; 1 2 3 4]
            };
    case 'infinite'
        data = {
            ShapeSet.LinearLine, ShapeSet.InfiniteLinearQuad, [1 2; 2 1]
            ShapeSet.LinearTria, ShapeSet.InfiniteLinearPenta, [1 2 3; 1 2 3]
            ShapeSet.LinearQuad, ShapeSet.InfiniteLinearHexa, [1 2 3 4; 1 2 3 4]
            };
end

Elem = [];

for i = 1 : size(data,1)
    srcLset = data{i,1};
    dstLset = data{i,2};
    idx = data{i,3};
    
    nDstNodes = size(dstLset.Nodes,1);
    
    sel = Elements(:,2) == srcLset.Id;
    if ~any(sel)
        continue;
    end
    
    elem = Elements(sel,:);
    nE = size(elem,1);
    
    dstElem = zeros(nRep*nE, 4+nDstNodes);
    for iRep = 1 : nRep
        dstElem((iRep-1)*nE+(1:nE),1:(4+nDstNodes)) =...
            [repmat([0 dstLset.Id],nE,1),...
            elem(:,3:4), ...
            elem(:,4+idx(1,:))+(iRep-1)*nNod, ...
            elem(:,4+idx(2,:))+iRep*nNod];
    end
    
    Elem(end+(1:size(dstElem,1)), 1:size(dstElem,2)) = dstElem;
end

end
