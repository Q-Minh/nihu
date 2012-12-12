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

%% Initialization
Elements = drop_IDs(model);
coords = model.Nodes(:,2:4);

cent = nan(size(Elements,1),3);
normal = nan(size(Elements,1),3);

for t = [12 23 24] % these element types are processed
    sel = Elements(:,2) == t;
    if any(sel)
        [cent(sel,:), normal(sel,:)] =...
            vert2gauss(1, coords, t, Elements(sel,4+(1:mod(t,10))));
    end
end

end
