function [cent, normal] = centnorm(model)
%CENTNORM Element centers and normals
%  [CENT, NORMAL] = CENTNORM(MESH) Returns the element centers CENT and the
%  element normals NORMAL of a surface mesh MESH.
%
% See also: vert2gauss, geo2gauss

%   Copyright 2008-2010 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modified: 2012.12.11.

%% Initialization
Elements = drop_IDs(model);
coords = model.Nodes(:,2:4);

cent = zeros(size(Elements,1),3);
normal = zeros(size(Elements,1),3);

for t = [12 23 24]
    sel = Elements(:,2) == t;
    if any(sel)
        [cent(sel,:), normal(sel,:)] =...
            vert2gauss(1, coords, t, Elements(sel,4+(1:mod(t,10))));
    end
end

end
