function Elements = flip_elements(Elements)
%FLIP_ELEMENTS  Flip elements of a NiHu mesh
%   MESH = FLIP_ELEMENTS(MESH) flips all the elements of a NiHu mesh

%   Copyright 2008-2010 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modified 2012.12.19.

%
LSetVector = [ShapeSet.LinearLine,...
    ShapeSet.LinearTria,...
    ShapeSet.LinearQuad, ...
    ShapeSet.LinearPenta, ...
    ShapeSet.LinearHexa
    ];

for i = 1 : length(LSetVector)
    lset = LSetVector(i);
    sel = find(Elements(:,2) == lset.Id);
    if ~isempty(sel)
        nNodes = size(lset.Nodes,1);
        ind = 4+(1:nNodes);
        Elements(sel,ind) = fliplr(Elements(sel,ind));
    end
end
