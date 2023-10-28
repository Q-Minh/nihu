function [GC, GN, W, GI, A] = geo2gauss(model, order)
%GEO2GAUSS Gaussian quadrature over a mesh
%   [XG, NG, WG, IG] = GEO2GAUSS(MESH, ORDER) computes  Gaussian integration
%   point coordinates XG, normal vectors NG and integration weights WG
%   over elements of a mesh.
% Parameters:
%   MESH : NiHu mesh structure.
%   ORDER : quadrature order.
%   XG : M x 3 matrix, Gaussian point coordinates
%   NG : M x 3 matrix, Unit normals in Gaussian points
%   WG : M x 1 vector, Gaussian integration weights
%   IG : M x 1 vector, element index of each Gaussian point
%   A  : matrix interpolating from nodes to Gaussian points
%
% Example:
%  sphere = create_sphere_boundary(1,10);
%  [xg, ng, wg] = geo2gauss(sphere, 1);
%  A = sum(w)
%
%  computes the area A of the unit sphere.
%
%  sphere = create_sphere_boundary(1,10);
%  [xg, ~, ~, ~, A] = geo2gauss(sphere, 3);
%  f = sin(2*pi*sphere.Nodes(:,2));
%  ox = A * f;
%
%  interpolates function 'x' to Gaussian coordinates
%
% See also: gaussquad1, gaussquad2, gaussquad3, centnorm, vert2gauss

%   Copyright 2008-2010 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modified: 2015.10.15.

% Initialization
Elements = drop_IDs(model);
coords = model.Nodes(:,2:4);

%
GC = zeros(0,3);
GN = zeros(0,3);
W = zeros(0,1);
GI = zeros(0,1);

% select unique lsets
lSets = ShapeSet.fromId(unique(Elements(:,2)));
I = [];
J = [];
V = [];
maxi = 0;
for iLset = 1 : length(lSets)
    lset = lSets(iLset);
    sel = Elements(:,2) == lset.Id;
    nNodes = size(lset.Nodes,1);
    [gc, gn, w, gi, an] = vert2gauss(order, coords, lset,...
        Elements(sel,4+(1:nNodes)));
    GC = [GC; gc]; %#ok<*AGROW>
    GN = [GN; gn];
    W = [W; w];
    GI = [GI; gi];
    if nargout > 4
        [i,j,v] = find(an);
        I = [I; i+maxi];
        maxi = max(I);
        J = [J; j];
        V = [V; v];
    end
end

% generate sparse interpolation matrix
if nargout > 4
    A = sparse(I, J, V, size(GC,1), size(coords,1));
end

end  %of function
