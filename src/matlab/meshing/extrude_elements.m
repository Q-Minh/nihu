function Elem = extrude_elements(mesh, nNod, nRep, reverse)
%EXTRUDE_ELEMENTS helper function for EXTRUDE_MESH and REVEOLVE_MESH
%   ELEM = EXTRUDE_ELEMENTS(MESH, NNOD, NREP, REVERSE)

% Create new elements
Elements = drop_IDs(mesh);

data = {
    ShapeSet.LinearLine, ShapeSet.LinearQuad, [1 2; 2 1]
    ShapeSet.LinearTria, ShapeSet.LinearPenta, [1 2 3; 1 2 3]
    ShapeSet.LinearQuad, ShapeSet.LinearHexa, [1 2 3 4; 1 2 3 4]
    };

Elem = zeros(nRep * size(Elements,1), 0);

for i = 1 : size(data,1)
    srcLset = data{i,1};
    dstLset = data{i,2};
    idx = data{i,3};
    
    nSrcNodes = size(srcLset.Nodes,1);
    nDstNodes = size(dstLset.Nodes,1);
    
    sel = Elements(:,2) == srcLset.Id;
    if ~any(sel)
        continue;
    end
    
    elem = Elements(sel,:);
    nE = size(elem,1);
    fl = reverse(sel);
    elem(fl,4+(1:nSrcNodes)) = fliplr(elem(fl,4+(1:nSrcNodes)));

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
