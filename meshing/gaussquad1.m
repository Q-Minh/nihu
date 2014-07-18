function [xx, ww] = gaussquad1(order)
%GAUSSQUAD1 1D Gaussian quadrature integration
%   [X, W] = GAUSSQUAD2(P) returns Gaussian integration base
%   points and weights for numerical integration over stadard 1D element.
% Parameters:
%   P     : Order of quadrature (2N-1) where N = number of Gauss points
%   X     : Gaussian quadrature base points
%   W     : Gaussian quadrature weights

%   Copyright 2008-2010 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modified 02.12.2009
%% Argument check
error(nargchk(1, 1, nargin, 'struct'));

%%
n = ceil((order+1)/2);
[xx, ww] = gaussquad(n); % Linear Gauss quad over (-1, +1)
end
