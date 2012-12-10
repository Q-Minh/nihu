function [cent, normal] = centnorm(model)
%CENTNORM Element centers and normals
%  [CENT, NORMAL] = CENTNORM(MESH) Returns the element centers CENT and the
%  element normals NORMAL of a surface mesh MESH.
%
% See also: vert2gauss, geo2gauss

%   Copyright 2008-2010 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modified: 18.11.2009.

%% Initialization
Elements = drop_IDs(model);
X = model.Nodes(:,2);
Y = model.Nodes(:,3);
Z = model.Nodes(:,4);

%% Process line elements
lin = find(Elements(:,2) == 12);
if ~isempty(lin)
    [cent(lin,:), normal(lin,:)] = vert2gauss(1, X,Y,Z, 12, Elements(lin,5:6));
end

%% Process tria elements
tri = find(Elements(:,2) == 23);
if ~isempty(tri)
    [cent(tri,:), normal(tri,:)] = vert2gauss(1, X,Y,Z, 23, Elements(tri,5:7));
end

%% Process quad elements
quad = find(Elements(:,2) == 24);
if ~isempty(quad)
    [cent(quad,:), normal(quad,:)] = vert2gauss(1, X,Y,Z, 24, Elements(quad,5:8));
end
end
