function x = burtmillgen(mesh, type, dist)
%BURTMILLGEN  Generate modified Burton and Miller points
%   X = BURTMILLGEN(MESH, DIST) generates modified Burton and Miller points
%   to numerically evaluate normal derivative of the BEM equations.
%
% Example:
%   s = create_brick_boundary(1,10);
%   x = burtmillgen(s);
%   plot3(x(:,1), x(:,2), x(:,3), '.');
%   plot_mesh(s);

% Last modified: 16.11.2009.
% Peter Fiala

if nargin < 2
    type = 'const';
end
% 1-point Gaussian quadrature;
[xg, ng, wg, ind] = geo2gauss(mesh, [1, 1]);
xg(ind,:) = xg;
ng(ind,:) = ng;
wg(ind) = wg;
if nargin < 3
    dist = sqrt(wg)/2;
end
x = xg + repmat(dist,1,3).*ng;