function model = flip_elements(model)
%FLIP_ELEMENTS  Flip elements of a NiHu mesh
%   MESH = FLIP_ELEMENTS(MESH) flips all the elements of a NiHu mesh

%   Copyright 2008-2010 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modified 2012.12.19.

%
types = [12 23 24 36 38];

for t = 1 : length(types)
    type = types(t);
    sel = find(model.Elements(:,2) == type);
    if ~isempty(sel)
        ind = 4+(1:mod(type, 10));
        model.Elements(sel,ind) = fliplr(model.Elements(sel,ind));
    end
end
