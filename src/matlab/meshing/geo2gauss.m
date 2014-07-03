function [gcoord, gnorm, weight, gind] = geo2gauss(model, order)
%GEO2GAUSS Gaussian quadrature over a mesh
%   [XG, NG, WG, IG] = GEO2GAUSS(MESH, P) computes  Gaussian integration
%   point coordinates XG, normal vectors NG and integration weights WG
%   over elements of a mesh.
% Parameters:
%   MESH : NiHu mesh structure.
%   P    : quadrature order.
%   XG : M x 3 matrix, Gaussian point coordinates
%   NG : M x 3 matrix, Unit normals in Gaussian points
%   WG : M x 1 vector, Gaussian integration weights
%   IG : M x 1 vector, element index of each Gaussian point
%
% Example:
%  sphere = create_sphere_boundary(1,10);
%  [xg, ng, wg] = geo2gauss(sphere, 1);
%  A = sum(w)
%
%  computes the area A of the unit sphere.
%
% See also: gaussquad1, gaussquad2, gaussquad3, centnorm, vert2gauss

%   Copyright 2008-2010 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modified: 2012.12.12.

%% Initialization
Elements = drop_IDs(model);
coords = model.Nodes(:,2:4);

%% Process line elements
lin = find(Elements(:,2) == 12);
if ~isempty(lin)
    [gc12, gn12, w12, gi12] = vert2gauss(order, coords, 12, Elements(lin,5:6));
    gi12 = lin(gi12);
else
    gc12 = zeros(0,3);
    gn12 = zeros(0,3);
    w12 = zeros(0,1);
    gi12 = zeros(0,1);
end

%% Process tria elements
tri = find(Elements(:,2) == 23);
if ~isempty(tri)
    [gc3, gn3, w3, gi3] = vert2gauss(order, coords, 23, Elements(tri,5:7));
    gi3 = tri(gi3);
else
    gc3 = zeros(0,3);
    gn3 = zeros(0,3);
    w3 = zeros(0,1);
    gi3 = zeros(0,1);
end

%% Process quad elements
quad = find(Elements(:,2) == 24);
if ~isempty(quad)
    [gc4, gn4, w4, gi4] = vert2gauss(order, coords, 24, Elements(quad,5:8));
    gi4 = quad(gi4);
else
    gc4 = zeros(0,3);
    gn4 = zeros(0,3);
    w4 = zeros(0,1);
    gi4 = zeros(0,1);
end

%% Combine TRIA and QUAD elements
gcoord = [gc12; gc3; gc4];
gnorm = [gn12; gn3; gn4];
weight = [w12; w3; w4];
gind = [gi12; gi3; gi4];
end
