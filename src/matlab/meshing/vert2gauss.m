function [gcoord, gnorm, weight, gind] = vert2gauss(order, coords, type, elem)
%VERT2GAUSS Gaussian quadrature from vertices and elements
%  [XG, NG, W, IG] = VERT2GAUSS(ORDER, COORDS, TYPE, ELEM)
%
% Example:
%
% See also: GEO2GAUSS, GAUSSQUAD1, GAUSSQUAD2, GAUSSQUAD3, SHAPEFUN

%   Copyright 2008-2010 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modified: 02.12.2009.

%%
nE = size(elem,1);  % number of elements

X = coords(:,1);
Y = coords(:,2);
Z = coords(:,3);

x = X(elem);        % vertices
y = Y(elem);
z = Z(elem);
% check if there was only one element
if size(x,2) == 1
    x = x.';
    y = y.';
    z = z.';
end
% Gaussian quadrature over the standard element
switch floor(type/10)
    case 1
        [xx, ww] = gaussquad1(order);
    case 2
        [xx, ww] = gaussquad2(order, size(elem,2));
    case 3
        [xx, ww] = gaussquad3(order, size(elem,2));
end
ng = length(ww);
[N, dN] = shapefun(xx, type);
% Gaussian quadrature coordinates
gx = N * x.';
gy = N * y.';
gz = N * z.';
gcoord = [gx(:) gy(:) gz(:)];
% Gaussian normals
switch floor(type/10)
    case 1
        gxx = dN * x.';
        gxy = dN * y.';
        gnorm = [gxy(:), -gxx(:), zeros(numel(gxx),1)];
    case 2
        gxx = dN(:,:,1) * x.';
        gxy = dN(:,:,1) * y.';
        gxz = dN(:,:,1) * z.';
        gyx = dN(:,:,2) * x.';
        gyy = dN(:,:,2) * y.';
        gyz = dN(:,:,2) * z.';
        gnorm = cross([gxx(:) gxy(:) gxz(:)], [gyx(:) gyy(:) gyz(:)], 2);
    case 3
        % TODO compute "normal" to get jac
        gnorm = zeros(size(gx,1), 3);
end
% Jacobian
j = sqrt(dot(gnorm,gnorm,2));
% Unit normal
gnorm = bsxfun(@times, gnorm, 1./j);
% Gaussian weight
weight = repmat(ww, nE, 1) .* j;
%
gind = repmat(1:nE,ng,1);
gind = gind(:);
end
