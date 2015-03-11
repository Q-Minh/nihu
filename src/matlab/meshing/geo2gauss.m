function [GC, GN, W, GI] = geo2gauss(model, order)
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

% Last modified: 2015.03.11.

%% Initialization
Elements = drop_IDs(model);
coords = model.Nodes(:,2:4);

%%
GC = zeros(0,3);
GN = zeros(0,3);
W = zeros(0,1);
GI = zeros(0,1);

for lset = [ShapeSet.LinearLine, ShapeSet.LinearTria, ShapeSet.LinearQuad]
    sel = Elements(:,2) == lset.Id;
    if any(sel)
        nNodes = size(lset.Nodes,1);
        [gc, gn, w, gi] = vert2gauss(order, coords, lset, Elements(sel,4+(1:nNodes)));
        GC = [GC; gc]; %#ok<*AGROW>
        GN = [GN; gn];
        W = [W; w];
        GI = [GI; gi];
    end
end

end
