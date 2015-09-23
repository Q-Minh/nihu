function Elem = extrude_elements(Elements, nNod, nRep, opt)
%EXTRUDE_ELEMENTS helper function for EXTRUDE_MESH and REVEOLVE_MESH
%   ELEM = EXTRUDE_ELEMENTS(MESH, NNOD, NREP, REVERSE)

% Create new elements

if nargin < 4
    opt = 'finite';
end

if nargin < 3
    nRep = 1;
end

switch opt
    case 'finite'
        extrudeRule = {
            ShapeSet.LinearLine, ShapeSet.LinearQuad, [1 2; 2 1]
            ShapeSet.LinearTria, ShapeSet.LinearPenta, [1 2 3; 1 2 3]
            ShapeSet.LinearQuad, ShapeSet.LinearHexa, [1 2 3 4; 1 2 3 4]
            };
    case 'infinite'
        extrudeRule = {
            ShapeSet.LinearLine, ShapeSet.InfiniteLinearQuad, [1 2; 2 1]
            ShapeSet.LinearTria, ShapeSet.InfiniteLinearPenta, [1 2 3; 1 2 3]
            ShapeSet.LinearQuad, ShapeSet.InfiniteLinearHexa, [1 2 3 4; 1 2 3 4]
            };
    otherwise
        error('NiHu:extrude_elements:Arg',...
            'invalid extrusion option ''%s''\nValid options are ''finite'' or ''infinite''',...
            opt);
end

Elem = [];

for i = 1 : size(extrudeRule,1)
    srcLset = extrudeRule{i,1};
    dstLset = extrudeRule{i,2};
    idx = extrudeRule{i,3};
    
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

end % of function
