function E = get_faces(elements)
%GET_FACES Extract all faces of a FE mesh
%   FACES = GET_FACES(ELEMENTS) extracts all the faces of a fe mesh given
%   by the matrix model.ELEMENTS.
%   The structure of the matrix ELEMENTS is
%   [ElemID LsetID MatID PropID NodID1 NODID2 ... NodIDn]
%   The structure of the matrix FACES is
%   [ElemID LsetID NodID1 ... NodIDn]
%   where ElemID refers to the parent element of the face.

%   Copyright 2008-2015 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

faceRules = {
    ShapeSet.LinearTria, ShapeSet.LinearLine, [1 2; 2 3; 3 1]
    ShapeSet.LinearQuad, ShapeSet.LinearLine, [1 2; 2 3; 3 4; 4 1]
    ShapeSet.LinearTetra, ShapeSet.LinearTria, [1 2 4; 2 3 4; 3 1 4; 1 3 2]
    ShapeSet.LinearPenta, ShapeSet.LinearTria, [1 2 3; 4 6 5]
    ShapeSet.LinearPenta, ShapeSet.LinearQuad, [1 2 5 4; 2 3 6 5; 3 1 4 6]
    ShapeSet.LinearHexa, ShapeSet.LinearQuad, [
    1 4 3 2; 5 6 7 8; 1 2 6 5;
    2 3 7 6; 3 4 8 7; 4 1 5 8 ]
    };

E = zeros(0,0);
N = 0;

for i = 1 : size(faceRules,1)
    srclset = faceRules{i,1};
    
    sel = elements(:,2) == srclset.Id;
    if ~any(sel)
        continue;
    end
    
    N = size(E,1);

    dstlset = faceRules{i,2};
    faceind = faceRules{i,3};
    
    f = reshape(faceind.', [], 1);
    e = elements(sel, 4+f);
    e = reshape(e.', size(faceind,2), []).';
    
    parentid = reshape(repmat(elements(sel,1).', size(faceind,1), 1), [], 1);
    
    ind = N + (1:size(parentid,1));
    
%     E(max(ind), 2+size(e,2)) = 0;   % zero padding
    E(ind,1) = parentid;
    E(ind,2) = dstlset.Id;
    E(ind,2+(1:size(e,2))) = e;
end

end
