function y = lintrans(x, X, Y)
%LINTRANS linear transformation from multi-interval to other
%   y = LINTRANS(x, X, Y) transforms x from multi-interval X to
%   multi-interval Y by means of a linear transform. Multi-intervals are ND
%   boxes, represented by a 2XN matrix. 
%
% Example:
%   y = lintrans([1 2; 1 3; 2 5], [0 0; 1 1], [0 0; 1 2]);

narginchk(3,3);

D = size(x,2);
if size(X,2) ~= D || size(Y,2) ~= D
    error('bbFMM:arg', 'bounding box and input dimensions must match');
end

% translate to origin
XC = bsxfun(@minus, x, mean(X,1));
% rescale
Z = bsxfun(@times, XC, diff(Y,1)./diff(X,1));
% translate to Y
y = bsxfun(@plus, Z, mean(Y, 1));

end
