function [cent, normal] = centnorm(model)
%CENTNORM Element centers and normals
%  [CENT, NORMAL] = CENTNORM(MESH) Returns the element centers CENT and the
%  element normals NORMAL of a line or surface mesh MESH. Only line, tria
%  and quad elements are processed. For other element types, nan vectors are
%  returned.
%
%  Example:
%    mesh = create_sphere_boundary(1,4);
%    [cent, normal] = centnorm(mesh);
%    hold on;
%    plot_mesh(mesh);
%    plot3(cent(:,1), cent(:,2), cent(:,3), 'r.');
%
% See also: VERT2GAUSS, GEO2GAUSS

%   Copyright 2008-2012 P. Fiala and P. Rucz
%   Budapest University of Technology and Economics

% Last modified: 2012.12.12.

% Initialization
Elements = drop_IDs(model);
coords = model.Nodes(:,2:4);

cent = nan(size(Elements,1),3);
normal = nan(size(Elements,1),3);

for lset = [ShapeSet.LinearLine ShapeSet.LinearTria ShapeSet.LinearQuad] % these element types are processed
    sel = Elements(:,2) == lset.Id;
    if any(sel)
        n = size(lset.Nodes,1);
        [cent(sel,:), normal(sel,:)] =...
            vert2gauss(1, coords, lset, Elements(sel,4+(1:n)));
    end
end

% Compute center for volume elements
for lset = [ShapeSet.LinearTetra ShapeSet.LinearPenta ShapeSet.LinearHexa]
    sel = Elements(:,2) == lset.Id;
    if any(sel)
        n = size(lset.Nodes,1);
        cent(sel,:) = vert2gauss(1, coords, lset, Elements(sel, 4+(1:n)));
        normal(sel,:) = zeros(numel(sel),3);
    end
end

end % of function
