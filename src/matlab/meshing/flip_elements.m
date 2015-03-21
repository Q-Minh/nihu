function Elements = flip_elements(Elements)
%FLIP_ELEMENTS  Flip elements of a NiHu mesh
%   MESH = FLIP_ELEMENTS(MESH) flips all the elements of a NiHu mesh

%   Copyright 2008-2010 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modified 2015.03.19.

% explains how a geometrical element is flipped
flipRules = {
    ShapeSet.ConstantLine.Id,  1
    ShapeSet.ConstantTria.Id,  1
    ShapeSet.ConstantQuad.Id,  1
    ShapeSet.ConstantPenta.Id,  1
    ShapeSet.ConstantHexa.Id,  1
    ShapeSet.LinearLine.Id,  [2 1]
    ShapeSet.LinearTria.Id,  [3 2 1]
    ShapeSet.LinearQuad.Id,  [4 3 2 1]
    ShapeSet.LinearTetra.Id, [3 2 1 4]
    ShapeSet.LinearPenta.Id,  [4 5 6 1 2 3]
    ShapeSet.LinearHexa.Id,  [5 6 7 8 1 2 3 4]
    ShapeSet.QuadraticLine.Id,  [3 2 1]
    ShapeSet.QuadraticTria.Id,  [3 2 1 6 5 4]
    };

% extract first column from flipRules
setids = cell2mat(flipRules(:,1));

% traverse unique lset id's
lsets = unique(Elements(:,2));
for i = 1 : length(lsets)
    id = lsets(i);
    lset = ShapeSet.fromId(id);
    nNodes = size(lset.Nodes,1);
    perm = flipRules{setids == id,2};
    sel = find(Elements(:,2) == id);
    ind = 4+(1:nNodes);
    nodes = Elements(sel,ind);
    Elements(sel,ind) = nodes(:,perm);
end

end % of function
