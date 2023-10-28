function [N, dN, ddN] = shapefun(xi, type)
%SHAPEFUN Shape functions over standard elements
% [N, dN] = SHAPEFUN(X, TYPE) returns linear shape function samples over
% standard elements. The input vector X is a nxd vector, each
% row containing the coordinates of a point in the standard
% element. The output N contains the n samples of the shape function
% N(xi), the output parameters dN and ddN contains the first and second
% derivatives for each location.

%   Copyright 2008-2015 P. Fiala and P. Rucz
%   Budapest University of Technology and Economics

% Last modified 2015.03.05.

s = ShapeSet.fromId(type);
[N, dN, ddN] = s.eval(xi);
end
