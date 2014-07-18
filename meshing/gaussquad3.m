function [xx, ww] = gaussquad3(order, nvert)
%GAUSSQUAD3 3D Gaussian quadrature integration
%   [X, W] = GAUSSQUAD3(N, Nvert) returns Gaussian integration base points
%   and weights for numerical integration over stadard 3D elements.
% Parameters:
%   N     : Number of Gaussian points.
%   Nvert : Number of vertices. 4 for standard TETRA, 6 for standard PENTA
%           and 8 for standard HEXA elements. The standard TETRA's
%           coordinates are
%           [(0,0,0), (1,0,0), (0,1,0), (0,0,1)]
%           The standard PENTA element's coordinates are
%           [(0,0,-1), (1,0,-1), (0,1,-1), (0,0,1), (1,0,1), (0,1,1)]
%           The standard HEXA element's coordinates are
%           [(-1,-1,1), (1,-1,1), (1,1,1), (-1,1,1),...
%            (-1,-1,1), (1,-1,1), (1,1,1), (-1,1,1)]
%   X     : Gaussian quadrature base points
%   W     : Gaussian quadrature weights
%
% See also: gaussquad, gaussquad2, vert2gauss, geo2gauss

%   Copyright 2008-2010 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modified 08.12.2009
%% Parameter check
error(nargchk(2, 2, nargin, 'struct'));

%%
switch nvert
    case 8
        [xi, w] = gaussquad1(order); % Linear Gauss quad over (-1, +1)
        [eta, xi, zeta] = meshgrid(xi,xi,xi);  % base points over unit cube
        W = (w * w.');                % weights over cube
        w = (W(:) * w.');             
        xx = [xi(:) eta(:) zeta(:)];
        ww = w(:);
    case 6
        [xi3, w3] = gaussquad2(order, 3);    % TRIA Gauss
        [xi2, w2] = gaussquad1(order);              % line Gauss
        x = repmat(xi3(:,1), 1, length(w2));        % merge the two
        y = repmat(xi3(:,2), 1, length(w2));
        z = repmat(xi2.', length(w3), 1);
        w = w3 * w2.';
        xx = [x(:) y(:) z(:)];
        ww = w(:);
    case 4
        error('NiHu:gaussquad3:TODO',...
            'Gaussian quadrature for TETRA not yet implemented!');
    otherwise
        error('NiHu:gaussquad3:argValue',...
            'Input argument Nvert should be 4 for TETRA, 6 for PENTA and 8 for HEXA');
end
end
